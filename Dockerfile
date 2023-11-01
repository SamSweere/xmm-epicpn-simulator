FROM ubuntu:latest

SHELL ["/bin/bash", "-c"]

ENV HOME=/xmm
WORKDIR $HOME
ENV SIMPUT=$HOME/simput SIXTE=$HOME/simput SAS_DIR=$HOME/sas SAS_CCFPATH=$HOME/ccf MINICONDA=$HOME/miniconda3

# Update everything and install required packages
RUN apt-get update && apt-get upgrade -y && apt-get dist-upgrade -y && \
    apt-get install software-properties-common -y &&  \
    apt-get update && apt-get upgrade -y && apt-get dist-upgrade -y && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && apt-get upgrade -y && apt-get dist-upgrade -y && \
    apt-get install -y git libtool autoconf wget rsync perl perlbrew libreadline-dev libncurses5-dev ncurses-dev curl \
    libcurl4 libcurl4-gnutls-dev xorg-dev make gcc g++ gfortran perl-modules libncurses-dev libexpat1-dev libgsl0-dev \
    libboost-dev libcmocka-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Download everything
RUN git clone http://www.sternwarte.uni-erlangen.de/git.public/simput.git/ && \
    git clone http://www.sternwarte.uni-erlangen.de/git.public/sixt/ && \
    wget https://www.sternwarte.uni-erlangen.de/~sixte/downloads/sixte/instruments/instruments_xmm-1.2.1.tar.gz && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    wget ftp://anonymous@sasdev-xmm.esac.esa.int/pub/sas/21.0.0/Linux/Ubuntu22.04/sas_21.0.0-Ubuntu22.04.tgz && \
    rsync -v -a --delete --delete-after --force --include='*.CCF' --exclude='*/' sasdev-xmm.esac.esa.int::XMM_VALID_CCF $SAS_CCFPATH && \
    wget https://heasarc.gsfc.nasa.gov/FTP/software/lheasoft/lheasoft6.32.1/heasoft-6.32.1src.tar.gz
    
# Instal SIMPUT and SIXTE
RUN cd $HOME/simput && autoreconf --install --force && \
    ./configure --prefix=$SIMPUT && make && make install && \
    cd $HOME/sixt && autoreconf --install --force && \
    ./configure --prefix=$SIXTE && make && make install && cd ~ && \
    # Add the SIXTE instrument files
    mv instruments_xmm-1.2.1.tar.gz $SIXTE &&  \
    cd $SIXTE && \
    tar zxf instruments_xmm-1.2.1.tar.gz &&  \
    rm instruments_xmm-1.2.1.tar.gz

# Install miniconda and initialize it
RUN bash miniconda.sh -b -u -p $MINICONDA && rm miniconda.sh
ENV PATH="$MINICONDA/bin:$PATH"
RUN conda init bash
RUN conda config --add channels conda-forge

# Create the conda environment
RUN conda create -n xmm python=3.11.5 astropy numpy matplotlib requests beautifultable scipy pypdf notebook astroquery \
    lxml yt h5py loguru

# Install perl-5.36.1 and the required modules
RUN mkdir -p $HOME/perl5/perlbrew/dists &&  \
    perlbrew --notest install perl-5.36.1 &&  \
    perlbrew switch perl-5.36.1 && \
    curl -L https://cpanmin.us | perl - App::cpanminus && \
    cpanm Switch && cpanm Shell && cpanm CGI

# Install XMM-SAS
ENV SAS_PERL=$HOME/perl5/perlbrew/perls/perl-5.36.1/bin/perl SAS_PYTHON=$MINICONDA/envs/xmm/bin/python
RUN mkdir $SAS_DIR && mv sas_21.0.0-Ubuntu22.04.tgz $SAS_DIR && cd $SAS_DIR && \
    tar zxf sas_21.0.0-Ubuntu22.04.tgz && \
    rm sas_21.0.0-Ubuntu22.04.tgz && \
    ./install.sh

# Install heasoft
ENV CC=/usr/bin/gcc CXX=/usr/bin/g++ FC=/usr/bin/gfortran PERL=$SAS_PERL PYTHON=$MINICONDA/envs/xmm/bin/python
RUN tar zxf heasoft-6.32.1src.tar.gz && \
    rm heasoft-6.32.1src.tar.gz && \
    unset CFLAGS CXXFLAGS FFLAGS LDFLAGS && \
    cd heasoft-6.32.1/BUILD_DIR/ && \
    ./configure && \
    make && \
    make install

# Set the environment variables for the image
RUN echo "export HOME=$HOME" >> $HOME/.bashrc && \
    echo "export SIMPUT=$SIMPUT" >> $HOME/.bashrc && \
    echo "export SIXTE=$SIXTE" >> $HOME/.bashrc && \
    echo "export PATH=\"$MINICONDA/bin:$HOME/perl5/perlbrew/perls/perl-5.36.1/bin/perl:\$PATH\"" >> $HOME/.bashrc && \
    echo "export SAS_DIR=$SAS_DIR/xmmsas_20230412_1735" >> $HOME/.bashrc && \
    echo "export SAS_CCF=$SAS_CCFPATH/ccf.cif" >> $HOME/.bashrc && \
    echo "export LD_LIBRARY_PATH=\"$MINICONDA/envs/xmm/lib:\$LD_LIBRARY_PATH\"" >> $HOME/.bashrc&& \
    echo "export HEADAS=`echo $HOME/heasoft-6.32.1/x*`" >> $HOME/.bashrc && \
    echo "export PYTHONPATH=\"\$SAS_DIR/lib/python:\$HEADAS/lib/python:\$PYTHONPATH\"" >> $HOME/.bashrc && \
    echo "export HEADASNOQUERY=" >> $HOME/.bashrc && \
    echo "export HEADASPROMPT=/dev/null" >> $HOME/.bashrc

# Start all installed software to make sure that everthing went fine and create ccf.cif
RUN . $HOME/.bashrc && \
    . $HEADAS/headas-init.sh && \
    . $SAS_DIR/setsas.sh && \
    . $SIXTE/bin/sixte-install.sh && \
    cd $SAS_CCFPATH && \
    cifbuild withobservationdate=yes
