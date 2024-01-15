FROM ubuntu:22.04

ENV HOME=/home/studtodorkov
RUN adduser studtodorkov --uid 1194 --disabled-password --gecos "" --home $HOME

ENV SIMPUT=$HOME/simput SAS_ROOT=$HOME/sas SAS_CCFPATH=$HOME/ccf MINICONDA=$HOME/miniconda3
ENV SAS_DIR=$SAS_ROOT/xmmsas_20230412_1735 SAS_CCF=$SAS_CCFPATH/ccf.cif HEADAS=$HOME/headas SIXTE=$SIMPUT

# Update everything and install required packages
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]
RUN apt-get update && apt-get upgrade -y && apt-get dist-upgrade -y && \
    apt-get install software-properties-common -y &&  \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && apt-get upgrade -y && apt-get dist-upgrade -y && \
    apt-get install -y git libtool autoconf wget rsync perl libreadline-dev libncurses5-dev ncurses-dev curl \
    libcurl4 libcurl4-gnutls-dev xorg-dev make gcc g++ gfortran perl-modules libncurses-dev libexpat1-dev libgsl0-dev \
    libboost-dev libcmocka-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

USER 1194
WORKDIR $HOME
COPY --chown=studtodorkov: --chmod=777 entrypoint.sh $HOME/
# Install miniconda and initialize it
COPY --chown=studtodorkov: --chmod=777 downloads/miniconda.sh $HOME/
RUN /bin/bash miniconda.sh -b -u -p ${MINICONDA} && rm miniconda.sh
ENV PATH="$MINICONDA/bin:$PATH"

RUN conda init bash
RUN conda create -c conda-forge -n xmm python=3.11.5 astropy numpy matplotlib requests beautifultable scipy  \
    pypdf notebook astroquery lxml yt h5py loguru pydantic jsonschema

# Install perl-5.36.1 and the required modules
RUN curl -L https://install.perlbrew.pl | bash
ENV PERLBREW_ROOT=$HOME/perl5/perlbrew
RUN $HOME/perl5/perlbrew/bin/perlbrew init && $HOME/perl5/perlbrew/bin/perlbrew install-cpanm
RUN mkdir -p ${HOME}/perl5/perlbrew/dists
RUN source ${HOME}/perl5/perlbrew/etc/bashrc && \
    $HOME/perl5/perlbrew/bin/perlbrew install -n -f perl-5.36.1
RUN $HOME/perl5/perlbrew/bin/perlbrew switch perl-5.36.1 && \
    $HOME/perl5/perlbrew/bin/cpanm -n Switch && $HOME/perl5/perlbrew/bin/cpanm -n Shell && $HOME/perl5/perlbrew/bin/cpanm -n CGI

COPY --chown=studtodorkov: --chmod=777 downloads/simput $HOME/simput_git/
RUN cd ${HOME}/simput_git && \
    echo "Initializing simput..." && autoreconf --install --force > /dev/null 2>&1 && \
    echo "Configuring simput..." && ./configure --prefix=${SIMPUT} > /dev/null 2>&1 &&  \
    echo "Building simput..." && make > /dev/null 2>&1 &&  \
    echo "Installing simput..." && make install > /dev/null 2>&1 && \
    echo "Cleaning simput..." && make clean > /dev/null 2>&1 && rm -rf ${HOME}/simput_git

COPY --chown=studtodorkov: --chmod=777 downloads/sixt $HOME/sixte_git/
RUN cd ${HOME}/sixte_git && \
    echo "Initializing sixte..." && autoreconf --install --force > /dev/null 2>&1 && \
    echo "Configuring sixte..." && ./configure --prefix=${SIXTE} > /dev/null 2>&1 &&  \
    echo "Building sixte..." && make > /dev/null 2>&1 &&  \
    echo "Installing sixte..." && make install > /dev/null 2>&1 && \
    echo "Cleaning sixte..." && make clean > /dev/null 2>&1 && rm -rf ${HOME}/sixte_git

WORKDIR $SIXTE
COPY --chown=studtodorkov: --chmod=777 downloads/instruments_xmm-1.2.1.tar.gz $SIXTE/
RUN tar zxf instruments_xmm-1.2.1.tar.gz && rm instruments_xmm-1.2.1.tar.gz

WORKDIR $SAS_ROOT
ENV SAS_PERL=$HOME/perl5/perlbrew/perls/perl-5.36.1/bin/perl SAS_PYTHON=$MINICONDA/envs/xmm/bin/python
COPY --chown=studtodorkov: --chmod=777 downloads/sas_21.0.0-Ubuntu22.04.tgz $SAS_ROOT/
USER 0
RUN tar zxf sas_21.0.0-Ubuntu22.04.tgz -C $SAS_ROOT && rm sas_21.0.0-Ubuntu22.04.tgz
RUN chown -R studtodorkov: $SAS_ROOT && chmod -R 777 $SAS_ROOT
USER 1194
RUN ./install.sh

COPY --chown=studtodorkov: --chmod=777 downloads/ccf $SAS_CCFPATH/

WORKDIR $HOME
ENV CC=/usr/bin/gcc CXX=/usr/bin/g++ FC=/usr/bin/gfortran PERL=$SAS_PERL PYTHON=$MINICONDA/envs/xmm/bin/python
COPY --chown=studtodorkov: --chmod=777 downloads/heasoft-6.32.1src.tar.gz $HOME/
RUN tar zxf heasoft-6.32.1src.tar.gz && rm heasoft-6.32.1src.tar.gz
RUN unset CFLAGS CXXFLAGS FFLAGS LDFLAGS && \
    cd heasoft-6.32.1/BUILD_DIR/ && \
    echo "Configuring heasoft..." && ./configure --prefix=${HEADAS} > /dev/null 2>&1 && \
    echo "Building heasoft..." && make > /dev/null 2>&1 && \
    echo "Installing heasoft..." && make install > /dev/null 2>&1 && \
    echo "Cleaning heasoft..." && make clean > /dev/null 2>&1 && \
    /bin/bash -c 'cd /home/studtodorkov/headas; for loop in x86_64*/*; do ln -sf $loop; done' \
    cd ${HOME}/heasoft-6.32.1 \
    cp -p Xspec/BUILD_DIR/hmakerc ${HEADAS}/bin/ \
    cp -p Xspec/BUILD_DIR/Makefile-std ${HEADAS}/bin/ \
    rm -rf Xspec/src/spectral


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

ENTRYPOINT ["/bin/bash", "${HOME}/entrypoint.sh"]
