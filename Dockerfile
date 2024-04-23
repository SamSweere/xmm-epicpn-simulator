FROM ubuntu:22.04 as builder

ENV HOME=/home/xmm_user
# Create user
RUN adduser xmm_user --uid 1000 --disabled-password --gecos "" --home $HOME

ENV SIMPUT=$HOME/simput SAS_ROOT=$HOME/sas SAS_CCFPATH=$HOME/ccf MINICONDA=$HOME/miniconda3
ENV SAS_DIR=$SAS_ROOT/xmmsas_20230412_1735 SAS_CCF=$SAS_CCFPATH/ccf.cif HEADAS=$HOME/headas SIXTE=$SIMPUT

# Update everything and install required packages
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]
RUN apt-get update && apt-get upgrade -y && apt-get dist-upgrade -y && \
    apt-get install software-properties-common -y &&  \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && apt-get upgrade -y && apt-get dist-upgrade -y && \
    apt-get install -y \
    # Heasoft packages
    build-essential curl gcc g++ gfortran libcurl4 libcurl4-gnutls-dev \
    libncurses5-dev libreadline6-dev \
    libcmocka-dev libexpat1-dev libgsl0-dev libfile-which-perl \
    libdevel-checklib-perl make ncurses-dev perl perl-modules xorg-dev \
    # SIMPUT
    autoconf libboost-dev libtool && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /var/cache/apt/*

USER 1000
WORKDIR $HOME

# Install miniconda and initialize it
COPY --chown=xmm_user: --chmod=777 downloads/miniconda.sh $HOME/
RUN /bin/bash miniconda.sh -b -u -p ${MINICONDA} && rm miniconda.sh
ENV PATH="$MINICONDA/bin:$PATH"

# Create the conda environment for xmm, also install numpy already since headas needs it
RUN conda init bash && \
    conda create -y -n xmm python=3.11.8 numpy astropy scipy matplotlib

# Build and install SIMPUT
COPY --chown=xmm_user: --chmod=777 downloads/simput_git $HOME/simput_git/
RUN cd ${HOME}/simput_git && \
    echo "Initializing simput..." && autoreconf --install --force && \
    echo "Configuring simput..." && ./configure --prefix=${SIMPUT} &&  \
    echo "Building simput..." && make &&  \
    echo "Installing simput..." && make install && \
    rm -rf ${HOME}/simput_git

# Build and install SIXTE
COPY --chown=xmm_user: --chmod=777 downloads/sixte_git $HOME/sixte_git/
RUN cd ${HOME}/sixte_git && \
    echo "Initializing sixte..." && autoreconf --install --force && \
    echo "Configuring sixte..." && ./configure --prefix=${SIXTE} &&  \
    echo "Building sixte..." && make &&  \
    echo "Installing sixte..." && make install && \
    rm -rf ${HOME}/sixte_git

# Extract and setup xmm instruments
WORKDIR $SIXTE
COPY --chown=xmm_user: --chmod=777 downloads/instruments_xmm-1.2.1.tar.gz $SIXTE/
RUN tar zxf instruments_xmm-1.2.1.tar.gz && rm instruments_xmm-1.2.1.tar.gz

# Extract and setup SAS
WORKDIR $SAS_ROOT
ENV SAS_PERL=/usr/bin/perl SAS_PYTHON=$MINICONDA/envs/xmm/bin/python
COPY --chown=xmm_user: --chmod=777 downloads/sas_21.0.0-Ubuntu22.04.tgz $SAS_ROOT/
USER 0
RUN tar zxf sas_21.0.0-Ubuntu22.04.tgz -C $SAS_ROOT && rm sas_21.0.0-Ubuntu22.04.tgz && \
    chown -R xmm_user: $SAS_ROOT && chmod -R 777 $SAS_ROOT
USER 1000
RUN ./install.sh

# Copy the Sas files
COPY --chown=xmm_user: --chmod=777 downloads/ccf $SAS_CCFPATH/

# Remove lib/libtinfo.so and lib/libtinfo.so.6 from the anaconda environment since it breaks bash
RUN rm $MINICONDA/envs/xmm/lib/libtinfo.so && rm $MINICONDA/envs/xmm/lib/libtinfo.so.6

# Install Heasoft
WORKDIR $HOME
ENV CC=/usr/bin/gcc CXX=/usr/bin/g++ FC=/usr/bin/gfortran PERL=$SAS_PERL PYTHON=$MINICONDA/envs/xmm/bin/python
COPY --chown=xmm_user: --chmod=777 downloads/heasoft-6.32.1src.tar.gz $HOME/
RUN tar zxf heasoft-6.32.1src.tar.gz && rm heasoft-6.32.1src.tar.gz

WORKDIR $HOME/heasoft-6.32.1/BUILD_DIR/
RUN unset CFLAGS CXXFLAGS FFLAGS LDFLAGS && \
    echo "Configuring heasoft..." && ./configure --prefix=${HEADAS} && \
    echo "Building heasoft..." && make && \
    echo "Installing heasoft..." && make install && \
    /bin/bash -c 'cd /home/xmm_user/headas; for loop in x86_64*/*; do ln -sf $loop; done' && \
    rm -rf ${HOME}/heasoft-6.32.1

WORKDIR $HOME

# Set the environment variables for the image
ENV HEADASPROMPT=/dev/null HEADASNOQUERY="" PATH=$MINICONDA/bin:$PERL:$PATH
ENV LD_LIBRARY_PATH=$MINICONDA/envs/xmm/lib:$LD_LIBRARY_PATH
ENV PYTHONPATH=$SAS_DIR/lib/python:$HEADAS/lib/python:$PYTHONPATH


# Start all installed software to make sure that everthing went fine and create ccf.cif
RUN . ${HEADAS}/headas-init.sh && \
    . ${SAS_DIR}/setsas.sh && \
    . ${SIXTE}/bin/sixte-install.sh && \
    cd ${SAS_CCFPATH} && \
    cifbuild withobservationdate=yes

# Create the entrypoint
RUN echo '#!/bin/bash' > ${HOME}/entrypoint.sh && \
    # echo 'echo "Arguments received: $@"' >> ${HOME}/entrypoint.sh && \
    echo '. $HEADAS/headas-init.sh' >> ${HOME}/entrypoint.sh && \
    echo '. $SAS_DIR/setsas.sh' >> ${HOME}/entrypoint.sh && \
    echo '. $SIXTE/bin/sixte-install.sh' >> ${HOME}/entrypoint.sh && \
    echo 'bash' >> ${HOME}/entrypoint.sh && \
    chmod +x ${HOME}/entrypoint.sh

# Add conda activate to the bashrc
RUN echo "source $MINICONDA/etc/profile.d/conda.sh && conda activate xmm" >> $HOME/.bashrc

# Stage 2: Final image
# We use the two stage build to keep the final image as small as possible (also to reduce layer caching)
FROM ubuntu:22.04

ENV HOME=/home/xmm_user

RUN adduser xmm_user --uid 1000 --disabled-password --gecos "" --home $HOME

# Set up environment variables
ENV HEADASPROMPT=/dev/null HEADASNOQUERY="" PATH=$MINICONDA/bin:$PERL:$PATH
ENV LD_LIBRARY_PATH=$MINICONDA/envs/xmm/lib:$LD_LIBRARY_PATH
ENV PYTHONPATH=$SAS_DIR/lib/python:$HEADAS/lib/python:$PYTHONPATH
ENV SIMPUT=$HOME/simput SAS_ROOT=$HOME/sas SAS_CCFPATH=$HOME/ccf MINICONDA=$HOME/miniconda3
ENV SAS_DIR=$SAS_ROOT/xmmsas_20230412_1735 SAS_CCF=$SAS_CCFPATH/ccf.cif HEADAS=$HOME/headas SIXTE=$SIMPUT

# Update install the runtime dependencies
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]
RUN apt-get update && apt-get upgrade -y && apt-get dist-upgrade -y && \
    apt-get install software-properties-common -y &&  \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && apt-get upgrade -y && apt-get dist-upgrade -y && \
    apt-get install -y \
    # Heasoft packages
    build-essential curl gcc g++ gfortran libcurl4 libcurl4-gnutls-dev \
    libncurses5-dev libreadline6-dev \
    libcmocka-dev libexpat1-dev libgsl0-dev libfile-which-perl \
    libdevel-checklib-perl make ncurses-dev perl perl-modules xorg-dev \
    # SIMPUT
    autoconf libboost-dev libtool \
    # General tools
    nano vim && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /var/cache/apt/*

# Copy the home directory
COPY --from=builder --chown=xmm_user: $HOME $HOME

# Finally Install the required python packages, we do this in the last step to avoid rebuilding the full image
# every time we change the requirements
WORKDIR $HOME
COPY --chown=xmm_user: --chmod=777 requirements.txt $HOME/
RUN $MINICONDA/envs/xmm/bin/pip install -r requirements.txt

# Set the working directory and entrypoint
WORKDIR $HOME
USER 1000
ENTRYPOINT ["/bin/bash", "entrypoint.sh"]
