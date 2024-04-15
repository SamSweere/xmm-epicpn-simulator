# XMM EPIC-pn Simulator
You can run this code either locally or via Docker. My main case for the latter is to run the code on a [Kubernetes](https://kubernetes.io) (aka K8s) cluster, which uses [CephFS](https://docs.ceph.com/en/nautilus/cephfs/). Therefore, my focus on optimization lies on utilising multiprocessing and handling of large amounts of small files.

The use cases do not differ much: The `Dockerfile` is a summary of steps to be taken to set everything up locally.

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

## Local setup
A word of warning: The following steps and the code run on Ubuntu 22.04. I haven't tested (nor do I intend to) other operating systems. If you're using Windows you'll probably run into problems, because of the logging library ([`Loguru`](https://github.com/Delgan/loguru)) I'm using. For more information see e. g. [this GitHub issue](https://github.com/Delgan/loguru/issues/1064).

The local setup is basically a step-by-step reproduction of the Dockerfile:
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
export MINICONDA=${HOME}/miniconda3

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
mkdir -p ${MINICONDA} && \
mkdir -p ${HOME}/simput_git/ && cd ${HOME}/simput_git/ && git clone http://www.sternwarte.uni-erlangen.de/git.public/simput.git/ . && \
mkdir -p ${HOME}/sixte_git/ && cd ${HOME}/sixte_git/ && git clone http://www.sternwarte.uni-erlangen.de/git.public/sixt/ . && \
cd ${HOME} && \
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ${MINICONDA}/miniconda.sh && \
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

## Configuring the code
The good thing: You'll need to fill out `config.json` only once! Every step relies on this configuration file and everything will be done accordingly. This file is divided into `environment`, `energy`, `download`, `simput` and `simulation`:

#### environment
This gives paths and some general information to the code:
- `working_dir`: Directory where the files will be saved to, while there are being worked on
- `output_dir`: If you're not running the code on a K8s cluster as I am, you can set this directory to the same value as `working_dir`. If `output_dir` and `working_dir` are not the same, the code will create a [tarball](https://manpages.ubuntu.com/manpages/focal/en/man1/tar.1.html) from the data in `working_dir` and move it to `output_dir`. This way the slow transfer speed of CephFS for small files is circumvented.
- `log_dir`: Directory where all the logs will be created. The code rotates the files every hour and keeps at most three files.
- `debug`: Switches of multiprocessing and gives more logs.
- `verbose`: Controls how much should be logged.
- `fail_on_error`: Sometimes errors may happen, which are not necessarily critical. For example: The download of one file from the Illustris Project could file, but everything else might work fine. If `fail_on_error` is `true`, then this will crash the program.
- `overwrite`: Controls if already existing files will be overwritten. Setting this to `false` could be useful if you want to avoid overwriting data that you have created previously. I would recommend to keep this to `true` and move any previously created files to some other directory.
- `consume_data`: Leave this at `false` if `working_dir == output_dir`! Otherwise files will be deleted to soon. This is mainly for my use case of running the code on K8s.

#### energy
Set the energy boundaries in `keV`:
- `emin`
- `emax`

#### download
- `num_processes`: How many processes should be run asynchrounosly. Recommended: As many CPUs as you have.
- `top_n`: How many cutouts to download for given simulations.
- `resolutions`: For every downloaded cutout, create images in these given resolutions
- `snapshots`: A dictionary of to-be-used snapshots with the corresponding redshift (see e.g. [TNG100-1](https://www.tng-project.org/data/downloads/TNG100-1/)). The IllustrisTNG project has snapshots 0 - 99.
- `simulations`: What simulations to consider with which width. Available are (with all of their sub-resolutions): `TNG50`, `TNG100`, and `TNG300`. The width is given as a tuple of `int` and `str`. If you don't want to use one simulation, then just delete it out of `config.json`.
- `modes`: There are two modes to create [`FITS`](https://heasarc.gsfc.nasa.gov/docs/heasarc/fits.html): projection and slice. The values given in the list are the axis for which the projection/slicing should be done. Both support the same values (`x`, `y`, `z`). If you want to use only one of the modes, then leave the list of the other empty.

#### simput
- `num_processes`: How many processes should be run asynchrounosly.
- `filter`: What XMM filter to use. Available: `thin`, `thick`, `med`. Only relevant for mode `bkg` (see below)
- `zoom_range`: From what range to randomly chose a zoom factor.
- `sigma_b_range`: The brightness sample range. This is based on the std of 50ks background. I.e. `sigma_b = 10` will result in a brightness of 10 times the background at 50ks.
- `offset_std`: The standard deviation of the normal distribution of the offset location around the bore-sight
- `num_img_sample`: How many simputs to create for previously downloaded files.
- `modes`: For what modes to create simputs. Available modes: `img`, `agn`, `bkg` (short for background). Set the value to `0` if none should be created. The mode `img` supports `-1`, which will create simputs for _all_ of the previously downloaded files. The mode `bkg` only supports a boolean value (or 0 and 1 accordingly).
- `instruments`: Only relevant for the mode `bkg`: For what instruments should a background simput be created. Available instruments: `epn`, `emos1`, `emos2`.

#### simulation
- `num_processes`: How many processes should be run asynchrounosly.
- `instrument_names`: What instruments should be simulated. Available instruments: `epn`, `emos1`, `emos2`.
- `filter`: What XMM filter to use. Available: `thin`, `thick`, `med`.
- `res_mults`: What resolution multiplication to simulate, e.g., 1x, 2x, 4x, etc.
- `max_exposure`: Max exposure to be simulated.
- `modes`: For what modes to run the instrument simulations. Available modes: `img`, `agn`, `bkg` (short for background). Set the value to `0` if none should be created. The modes `img` and `agn` support `-1`, which will run the simulation for _all_ of the previously created simputs for that mode.
- `sim_separate_ccds`: If the individual CCDs of XMM should be simulated or if they should be considered as "one big CCD".
- `wait_time`: If not 0, then Out-Of-Time events will be simulated.

## Running the code
The code is split up into different steps, represented by different scripts. If you want to go through the whole process, then you _must_ execute the steps in the correct order. They are numbered accordingly. There are following steps:

1. `01_download_files.py`: Download files from the [Illustris Project](https://www.tng-project.org). Before you can do that you'll need an API key. For this check out their [registration page](https://www.tng-project.org/users/register/). After your request has been approved, you'll see your personal API key after you login. Please keep this key to yourself!

2. `02_generate_simput.py`: Create SIMPUT files based on the previously downloaded files.

3. `03_xmm_simulation.py`: Simulate XMM-Newton for the previously created SIMPUT files. TBD: I will add at least one other satellite to choose from.

4. `04_combine_simulations.py`: **Not used right now!** I will rewrite this step to merge images from different satellites/different sensors.

Executing any of the scripts is same for both setups:

1. Set your configuration parameters as needed (see above)
2. Initialise external tools:

```shell
. ${HEADAS}/headas-init.sh && . ${SAS_DIR}/setsas.sh && . ${SIXTE}/bin/sixte-install.sh
```

3. Choose what step you want to run
4. Run `conda run -n xmm --no-capture-output python /path/to/script` with the needed command line arguments:

   1. `01_download_files.py` requires two arguments:

       1. `-k` followed by your personal Illustris API key (see below)
       2. `-p` followed by the path to the `config.json`
   2. `02_generate_simput.py` requires three arguments:
       1. `-a` followed by the path to the `agn_counts.cgi` file in `res`
       2. `-p` followed by the path to the `config.json`
       3. `-s` followed by the path to `res/spectrums`
    3. `03_xmm_simulation.py` requires one argument:
       1. `-p` followed by the path to the `config.json`

## IllustrisTNG simulations
For our XMM simulations we need sources to simulate (simulation input).
In our project we are especially interested in extended sources.
We take these extended sources from the Illustris TNG project (https://www.tng-project.org/).
This is a large cosmological hydrodynamical simulation of galaxy formation containing hundreds
of terabytes of simulated data. From this we take the most massive objects and take
x-ray projections and x-ray slices (less realistic but contains more clearly defined structure).
Note that cutout files are relatively large (100-1000 mb) and can take a while to download, it will first download all the relevant cutouts before generating the images.

#### Notes on extended sources (images as simput)
To create the simput for extended sources we use fits image files. In order to have a realistic distribution we augment these images using:
- Brightness: The brightness of the source is internally defined as sigma_b. This is based on the std of 50ks background. I.e. `sigma_b = 10` will result in a brightness of 10 times the background at 50ks.
The images are used as a distribution of a given brightness.
We determine the final brightness by taking a center cutout of the image and set this to the brightness defined by sigma_b.
- Location: We augment to location by offsetting the image from the bore-axis. Since real xmm observation are usually focussed on the center of extended sources we by default offcenter the images by a small amount around the bore-sight based on a normal distribution.
- Size (zoom): We augment the size of the extended source by artificially zooming in or out.

## XMM Simulation
The XMM simulations are done using SIXTE X-ray simulation software (https://www.sternwarte.uni-erlangen.de/research/sixte/).
All the elements that make up a XMM observation are simulated separately: extended source, agn and background.
These can then in the future be combined with a detector-mask to create a realistic XMM observation.
Since this is a simulation we can also simulate observations where XMM has a higher resolution (both spatial and psf wise).

## Acknowledgements
