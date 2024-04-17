# Installation
This code can be executed either locally or through Docker. The primary use case for Docker is to deploy the code on a [Kubernetes](https://kubernetes.io) (K8s) cluster utilizing [CephFS](https://docs.ceph.com/en/nautilus/cephfs/). Optimization efforts are concentrated on implementing multiprocessing and efficiently managing large volumes of small files.

The deployment process is similar for both use cases: the `Dockerfile` provides a detailed set of instructions for setting up the environment, whether locally or on a cluster.


## Downloading Pre-built Docker Image
The easiest way to get started is to download the pre-built Docker image from the [GitHub Container Registry](https://hub.docker.com/repository/docker/samsweere/xmm-epicpn-simulator/general)
```shell
# TODO
```

Note that there are no guarantees that the image is up-to-date or available in the future. If the image is not available, you can build the docker image yourself, see the [Docker Setup](#docker-setup) section below.

## Docker Setup
Docker is a containerised system which makes it possible create and run a whole operation system in a container.
If not already done, you need to [install](https://docs.docker.com/engine/install/) the Docker engine. If you want a GUI, you can instead install [Docker Desktop](https://docs.docker.com/desktop/).

1. If not already done, pull this repo, open a terminal and switch directory into the root directory of the repo.
2. Download everything beforehand to reduce the image build time:

```shell
mkdir downloads/ && cd downloads/ && \
mkdir simput_git/ && cd simput_git/ && git clone http://www.sternwarte.uni-erlangen.de/git.public/simput.git/ . && cd ../ && \
mkdir sixte_git/ && cd sixte_git/ && git clone http://www.sternwarte.uni-erlangen.de/git.public/sixt/ . && cd ../ && \
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
wget https://www.sternwarte.uni-erlangen.de/~sixte/downloads/sixte/instruments/instruments_xmm-1.2.1.tar.gz && \
wget ftp://anonymous@sasdev-xmm.esac.esa.int/pub/sas/21.0.0/Linux/Ubuntu22.04/sas_21.0.0-Ubuntu22.04.tgz && \
rsync -v -a --delete --delete-after --force --include='*.CCF' --exclude='*/' sasdev-xmm.esac.esa.int::XMM_VALID_CCF ccf/ && \
wget https://heasarc.gsfc.nasa.gov/FTP/software/lheasoft/lheasoft6.32.1/heasoft-6.32.1src.tar.gz
```

3. Build the image given by the `Dockerfile` in this repo: `docker build -t {your_own_tag} .` This will take an hour or two.
4. If needed, login into the image registries used by your cluster via `docker login {image_registry_url}` and push the image via `docker push {your_own_tag}`
5. You can access the image via `docker run -it {your_own_tag}` and run the code as described below.

## Local Installation
A word of warning: The following steps and the code run on Ubuntu 22.04. I haven't tested (nor do I intend to) other operating systems. If you're using Windows you'll probably run into problems, because of the logging library ([`Loguru`](https://github.com/Delgan/loguru)) I'm using. For more information see e. g. [this GitHub issue](https://github.com/Delgan/loguru/issues/1064).

The local installation is basically a step-by-step reproduction of the Dockerfile:
### 1. Setup your Ubuntu:
```shell
sudo apt-get install software-properties-common && \
sudo add-apt-repository ppa:ubuntu-toolchain-r/test && \
sudo apt-get update && sudo apt-get upgrade -y && sudo apt-get dist-upgrade -y && \
sudo apt-get install -y git libtool autoconf wget rsync perl libreadline-dev \
libncurses5-dev ncurses-dev curl libcurl4 libcurl4-gnutls-dev xorg-dev make gcc \
g++ gfortran perl-modules libncurses-dev libexpat1-dev libgsl0-dev libboost-dev \
libcmocka-dev
```
Setup environment variables _permanentely_ for both login and non-login shells. Add following lines to `${HOME}/.profile`:

```shell
# You can adjust these values as you like. These are my preferred paths.
export SIMPUT=${HOME}/simput
export SIXTE=${SIMPUT} # This path has to be the same as SIMPUT!
# export MINICONDA=${HOME}/miniconda3

# SAS variables
export SAS_ROOT=${HOME}/sas
export SAS_DIR=${SAS_ROOT}/xmmsas_20230412_1735 # xmmsas_20230412_1735 is mandatory!
export SAS_CCFPATH=${HOME}/ccf
export SAS_CCF=${SAS_CCFPATH}/ccf.cif
export SAS_PERL=${HOME}/perl5/perlbrew/perls/perl-5.36.1/bin/perl
export SAS_PYTHON=${MINICONDA}/envs/xmm/bin/python

# HEADAS variables
export HEADAS=${HOME}/headas
export HEADASPROMPT=/dev/null
export HEADASNOQUERY=""
export LD_LIBRARY_PATH=${MINICONDA}/envs/xmm/lib:${LD_LIBRARY_PATH}

# General
export PATH=${MINICONDA}/bin:${PERL}:${PATH}
export PYTHONPATH=${SAS_DIR}/lib/python:${HEADAS}/lib/python:${PYTHONPATH}
```
Make sure that these variables are loaded. If you want to make sure you can restart your PC.

### 2. Download everything
This will take some time...
```shell
# mkdir -p ${MINICONDA} && \
mkdir -p ${HOME}/simput_git/ && cd ${HOME}/simput_git/ && git clone http://www.sternwarte.uni-erlangen.de/git.public/simput.git/ . && \
mkdir -p ${HOME}/sixte_git/ && cd ${HOME}/sixte_git/ && git clone http://www.sternwarte.uni-erlangen.de/git.public/sixt/ . && \
cd ${HOME} && \
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ${MINICONDA}/miniconda.sh && \
wget https://www.sternwarte.uni-erlangen.de/~sixte/downloads/sixte/instruments/instruments_xmm-1.2.1.tar.gz && \
wget ftp://anonymous@sasdev-xmm.esac.esa.int/pub/sas/21.0.0/Linux/Ubuntu22.04/sas_21.0.0-Ubuntu22.04.tgz && \
rsync -v -a --delete --delete-after --force --include='*.CCF' --exclude='*/' sasdev-xmm.esac.esa.int::XMM_VALID_CCF $SAS_CCFPATH && \
wget https://heasarc.gsfc.nasa.gov/FTP/software/lheasoft/lheasoft6.32.1/heasoft-6.32.1src.tar.gz
```


### 3. Setup python environment
There are two ways to setup the python environment: Either via `pyenv` or `miniconda`. I would recommend using `pyenv` and `poetry` since the dependencies are already defined in `pyproject.toml`. If you want to use `miniconda`, then you'll need to install the dependencies manually.

### 3. Setup `python` using `pyenv` and `poetry`
The following steps are for a detailed installation of the development environment. Note that for every step there are multiple ways to i.e. install python, create an environment or install dependencies. The following steps are just one way to do it.


1. Install `Pyenv`:
    https://github.com/pyenv/pyenv#installation
2. Install `python 3.11.8`:
    ```bash
    pyenv install 3.11.8
    ```
3. Install pyenv-virtualenv:
    https://github.com/pyenv/pyenv-virtualenv

4. Create a virtual environment:
    ```bash
    pyenv virtualenv 3.11.8 xmm-epicpn-simulator
    ```
5. Enable and use the virtual environment:
    ```bash
    pyenv local xmm-epicpn-simulator
    pyenv activate xmm-epicpn-simulator
    ```
6. Install poetry:
    ```bash
    pip install poetry
    ```
7. Install the dependencies:
    ```bash
    poetry install
    ```

### (Optionally if you want to use conda) 3. Setup `miniconda`
Install `miniconda`:
```shell
bash ${MINICONDA}/miniconda.sh -b -u -p ${MINICONDA} && rm -rf ${MINICONDA}/miniconda.sh
```
Initialise `conda`:
```shell
conda init bash
```
_Restart your terminal._

Create an conda environment:
```shell
conda create -c conda-forge -n xmm python=3.11.5 astropy numpy matplotlib requests beautifultable scipy pypdf notebook astroquery lxml yt h5py loguru pydantic
```
Activate our new conda environment:
```shell
conda activate xmm
```

### 4. Setup Perl
Install [`perlbrew`](https://perlbrew.pl):
```shell
\curl -L https://install.perlbrew.pl | bash
```
Initialise `perlbrew` and install `cpanm`:
```shell
perlbrew init && perlbrew install-cpanm
```
_Restart your terminal._

Install `perl-5.36.1` and needed modules
```shell
perlbrew install --switch -n -f perl-5.36.1 && cpanm -n Switch && cpanm -n Shell && cpanm -n CGI
```

### 5. Setup SIMPUT and SIXTE
Install `SIMPUT`:
```shell
cd ${HOME}/simput_git/ && \
echo "Initializing simput..." && autoreconf --install --force > /dev/null 2>&1 && \
echo "Configuring simput..." && ./configure --prefix=${SIMPUT} > /dev/null 2>&1 &&  \
echo "Building simput..." && make > /dev/null 2>&1 &&  \
echo "Installing simput..." && make install > /dev/null 2>&1 && \
echo "Cleaning simput..." && make clean > /dev/null 2>&1 && rm -rf ${HOME}/simput_git/ && \
cd ${HOME}
```
Install `SIXTE`:
```shell
cd ${HOME}/sixte_git/ && \
echo "Initializing sixte..." && autoreconf --install --force > /dev/null 2>&1 && \
echo "Configuring sixte..." && ./configure --prefix=${SIMPUT} > /dev/null 2>&1 &&  \
echo "Building sixte..." && make > /dev/null 2>&1 &&  \
echo "Installing sixte..." && make install > /dev/null 2>&1 && \
echo "Cleaning sixte..." && make clean > /dev/null 2>&1 && rm -rf ${HOME}/sixte_git/ && \
cd ${HOME}
```
Get the needed instrument files:
```shell
mv instruments_xmm-1.2.1.tar.gz ${SIXTE}/ && cd ${SIXTE} && \
tar zxf instruments_xmm-1.2.1.tar.gz && rm instruments_xmm-1.2.1.tar.gz && \
cd ${HOME}
```
### 6. Setup `SAS`
Install `SAS`:
```shell
mkdir -p ${SAS_ROOT} && mv sas_21.0.0-Ubuntu22.04.tgz ${SAS_ROOT}/ && cd ${SAS_ROOT} && \
tar zxf sas_21.0.0-Ubuntu22.04.tgz -C $SAS_ROOT && rm sas_21.0.0-Ubuntu22.04.tgz && ./install.sh && \
cd ${HOME}
```
### 7. Setup `headas`
```shell
tar zxf heasoft-6.32.1src.tar.gz && rm heasoft-6.32.1src.tar.gz && \
unset CFLAGS CXXFLAGS FFLAGS LDFLAGS && \
CC=/usr/bin/gcc CXX=/usr/bin/g++ FC=/usr/bin/gfortran PERL=${SAS_PERL} PYTHON=${MINICONDA}/envs/xmm/bin/python && \
cd heasoft-6.32.1/BUILD_DIR/ && \
echo "Configuring heasoft..." && ./configure --prefix=${HEADAS} > /dev/null 2>&1 && \
echo "Building heasoft..." && make > /dev/null 2>&1 && \
echo "Installing heasoft..." && make install > /dev/null 2>&1 && \
echo "Cleaning heasoft..." && make clean > /dev/null 2>&1 && \
/bin/bash -c 'cd ${HEADAS}; for loop in x86_64*/*; do ln -sf $loop; done' && \
cd ${HOME}/heasoft-6.32.1 && \
cp -p Xspec/BUILD_DIR/hmakerc ${HEADAS}/bin/ && \
cp -p Xspec/BUILD_DIR/Makefile-std ${HEADAS}/bin/ && \
rm -rf Xspec/src/spectral && \
cd ${HOME}
```

### 8. Check if everything works fine
_Restart your terminal_ and run:
```shell
. ${HEADAS}/headas-init.sh && \
. ${SAS_DIR}/setsas.sh && \
. ${SIXTE}/bin/sixte-install.sh && \
cd ${SAS_CCFPATH} && \
cifbuild withobservationdate=yes && \
cd ${HOME}
```
If that runs fine, you're ready to go!
