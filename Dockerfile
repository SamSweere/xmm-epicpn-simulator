FROM ubuntu:latest

SHELL ["/bin/bash", "-c"]

ENV HOME=/xmm
WORKDIR $HOME

# Update everything and install required packages
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y git libtool autoconf wget rsync perl perlbrew && \
    apt-get install -y libreadline-dev libncurses5-dev ncurses-dev curl libcurl4 libcurl4-gnutls-dev xorg-dev make  \
    gcc g++ gfortran perl-modules && \
    apt-get install -y libncurses-dev libexpat1-dev libgsl0-dev libboost-dev libcmocka-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Instal SIMPUT and SIXTE from source
RUN git clone http://www.sternwarte.uni-erlangen.de/git.public/simput.git/ && \
    cd simput && autoreconf --install --force && \
    ./configure --prefix=/xmm/simput && make && make install && \
    git clone http://www.sternwarte.uni-erlangen.de/git.public/sixt/ && \
    cd sixt && autoreconf --install --force && \
    ./configure --prefix=/xmm/simput && make && make install && \
    echo "export SIMPUT=$HOME/simput" >> ~/.bashrc && \
    echo "export SIXTE=$HOME/simput" >> ~/.bashrc

ENV SIMPUT=$HOME/simput SIXTE=$HOME/simput

# Check that the installation of SIXTE was successfull
RUN . $SIXTE/bin/sixte-install.sh

# Download the SIXTE instrument files
RUN cd $SIXTE && \
    wget https://www.sternwarte.uni-erlangen.de/~sixte/downloads/sixte/instruments/instruments_xmm-1.2.1.tar.gz && \
    tar zxf instruments_xmm-1.2.1.tar.gz && rm instruments_xmm-1.2.1.tar.gz

# Install miniconda and initialize it
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -u -p $HOME/miniconda3 && \
    rm miniconda.sh && \
    echo "export PATH=$HOME/miniconda3/bin:$PATH" >> ~/.bashrc
ENV PATH="$HOME/miniconda3/bin:$PATH"
RUN conda init bash
RUN conda config --add channels conda-forge

# Create the conda environment
RUN conda create -n xmm python=3.11.5 astropy numpy matplotlib requests beautifultable scipy pypdf notebook astroquery tqdm lxml

ENV SAS_PYTHON=$HOME/miniconda3/envs/xmm/bin/python

# Install perl-5.36.1 and the required modules
RUN mkdir -p ~/perl5/perlbrew/dists &&  \
    perlbrew --notest install perl-5.36.1 &&  \
    perlbrew switch perl-5.36.1 && \
    cpan App::cpanminus && cpan Module::Switch && cpan Module::Shell && cpan Module::CGI

ENV SAS_PERL=~/perl5/perlbrew/perls/perl-5.36.1/bin/perl

# Install XMM-SAS
RUN mkdir $HOME/sas && cd $HOME/sas && \
    wget ftp://anonymous@sasdev-xmm.esac.esa.int/pub/sas/21.0.0/Linux/Ubuntu22.04/sas_21.0.0-Ubuntu22.04.tgz && \
    tar zxf sas_21.0.0-Ubuntu22.04.tgz && \
    rm sas_21.0.0-Ubuntu22.04.tgz && \
    ./install.sh && \
    echo "SAS_DIR=$HOME/sas/xmmsas_20230412_1735" >> ~/.bashrc

ENV SAS_DIR=$HOME/sas/xmmsas_20230412_1735

# Get the CCF and set SAS_CCFPATH
RUN mkdir $HOME/xmm_ccf && cd $HOME/xmm_ccf &&  \
    rsync -v -a --delete --delete-after --force --include='*.CCF' --exclude='*/' sasdev-xmm.esac.esa.int::XMM_VALID_CCF . && \
    echo "SAS_CCFPATH=$HOME/xmm_ccf" >> ~/.bashrc

ENV SAS_CCFPATH=$HOME/xmm_ccf

# Install heasoft
ENV CC=/usr/bin/gcc CXX=/usr/bin/g++ FC=/usr/bin/gfortran PERL=~/perl5/perlbrew/perls/perl-5.36.1/bin/perl PYTHON=$HOME/miniconda3/envs/xmm/bin/python
RUN wget https://heasarc.gsfc.nasa.gov/FTP/software/lheasoft/lheasoft6.32.1/heasoft-6.32.1src.tar.gz && \
    tar zxf heasoft-6.32.1src.tar.gz && \
    rm heasoft-6.32.1src.tar.gz && \
    unset CFLAGS CXXFLAGS FFLAGS LDFLAGS && \
    cd heasoft-6.32.1/BUILD_DIR/ && \
    ./configure && \
    make && \
    make install && \
    echo "export HEADAS=`echo $HOME/heasoft-6.32.1/x*`" >> ~/.bashrc

ENV PYTHONPATH="$PYTHONPATH:$SAS_DIR/lib/python:$HEADAS/lib/python"

# Check if HEASoft, XMM-SAS have been installed correctly
RUN . $HEADAS/headas-init.sh && \
   . $SAS_DIR/setsas.sh
