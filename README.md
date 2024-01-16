# XMM EPIC-pn Simulator
You can run this code either locally or via Docker. My main case for the latter is to run the code on a [Kubernetes](https://kubernetes.io) (aka K8s) cluster, which uses [CephFS](https://docs.ceph.com/en/nautilus/cephfs/). Therefore, my focus on optimization lies on utilising multiprocessing and handling of large amounts of small files.

The use cases do not differ much: The `Dockerfile` is a summary of steps to be taken to set everything up locally. 

## Docker Setup
Docker is a containerised system which makes it possible create and run a whole operation system in a container.
If not already done, you need to [install](https://docs.docker.com/engine/install/) the Docker engine. If you want a GUI, you can instead install [Docker Desktop](https://docs.docker.com/desktop/).

_Download everything beforehand_

1. If not already done, pull this repo, open a terminal and switch directory into the root directory of the repo.
2. Build the image given by the `Dockerfile` in this repo: `docker build -t {your_own_tag} .` This will take an hour or two.
3. If needed, login into the image registries used by your cluster via `docker login {image_registry_url}` and push the image via `docker push {your_own_tag}`

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
Setup environment variables _permanentely_. Please research how you can do that for your shell (see e. g. [here](https://unix.stackexchange.com/questions/117467/how-to-permanently-set-environmental-variables)). If you're using `bash` add following lines to e. g. `${HOME}/.bashrc`:

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
mkdir -p ${MINICONDA} \
mkdir -p ${HOME}/simput_git/ && cd ${HOME}/simput_git/ && git clone http://www.sternwarte.uni-erlangen.de/git.public/simput.git/ . \
mkdir -p ${HOME}/sixte_git/ && cd ${HOME}/sixte_git/ && git clone http://www.sternwarte.uni-erlangen.de/git.public/sixt/ . \
cd ${HOME} \
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ${MINICONDA}/miniconda.sh \
wget https://www.sternwarte.uni-erlangen.de/~sixte/downloads/sixte/instruments/instruments_xmm-1.2.1.tar.gz \
wget ftp://anonymous@sasdev-xmm.esac.esa.int/pub/sas/21.0.0/Linux/Ubuntu22.04/sas_21.0.0-Ubuntu22.04.tgz \
rsync -v -a --delete --delete-after --force --include='*.CCF' --exclude='*/' sasdev-xmm.esac.esa.int::XMM_VALID_CCF $SAS_CCFPATH \
wget https://heasarc.gsfc.nasa.gov/FTP/software/lheasoft/lheasoft6.32.1/heasoft-6.32.1src.tar.gz
```
### 3. Setup `miniconda`
Install `miniconda`:
```shell
. ${MINICONDA}/miniconda.sh -b -u -p ${MINICONDA} && rm -rf ${MINICONDA}/miniconda.sh
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
cd ${HOME}/simput_git/ \
echo "Initializing simput..." && autoreconf --install --force > /dev/null 2>&1 && \
echo "Configuring simput..." && ./configure --prefix=${SIMPUT} > /dev/null 2>&1 &&  \
echo "Building simput..." && make > /dev/null 2>&1 &&  \
echo "Installing simput..." && make install > /dev/null 2>&1 && \
echo "Cleaning simput..." && make clean > /dev/null 2>&1 && rm -rf ${HOME}/simput_git/ \
cd ${HOME}
```
Install `SIXTE`:
```shell
cd ${HOME}/sixte_git/ \
echo "Initializing sixte..." && autoreconf --install --force > /dev/null 2>&1 && \
echo "Configuring sixte..." && ./configure --prefix=${SIMPUT} > /dev/null 2>&1 &&  \
echo "Building sixte..." && make > /dev/null 2>&1 &&  \
echo "Installing sixte..." && make install > /dev/null 2>&1 && \
echo "Cleaning sixte..." && make clean > /dev/null 2>&1 && rm -rf ${HOME}/sixte_git/ \
cd ${HOME}
```
Get the needed instrument files:
```shell
mv instruments_xmm-1.2.1.tar.gz ${SIXTE}/ && cd ${SIXTE} \
tar zxf instruments_xmm-1.2.1.tar.gz && rm instruments_xmm-1.2.1.tar.gz \
cd ${HOME}
```
### 6. Setup `SAS`
Install `SAS`:
```shell
mkdir -p ${SAS_ROOT} && mv sas_21.0.0-Ubuntu22.04.tgz ${SAS_ROOT}/ && cd ${SAS_ROOT} \
tar zxf sas_21.0.0-Ubuntu22.04.tgz -C $SAS_ROOT && rm sas_21.0.0-Ubuntu22.04.tgz && ./install.sh \
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
/bin/bash -c 'cd ${HEADAS}; for loop in x86_64*/*; do ln -sf $loop; done' \
cd ${HOME}/heasoft-6.32.1 \
cp -p Xspec/BUILD_DIR/hmakerc ${HEADAS}/bin/ \
cp -p Xspec/BUILD_DIR/Makefile-std ${HEADAS}/bin/ \
rm -rf Xspec/src/spectral
```

### 8. Check if everything works fine
_Restart your terminal_ and run:
```shell
. ${HEADAS}/headas-init.sh && \
. ${SAS_DIR}/setsas.sh && \
. ${SIXTE}/bin/sixte-install.sh && \
cd ${SAS_CCFPATH} && \
cifbuild withobservationdate=yes
```
If that runs fine, you're ready to go!

## Running the code
The code is split up into different steps, represented by different scripts. If you want to go through the whole process, then you _must_ execute the steps in the correct order. They are numbered accordingly. There are following steps:

1. `01_download_files.py`: Download files from the [Illustris Project](https://www.tng-project.org). Before you can do that you'll need an API key. For this check out their [registration page](https://www.tng-project.org/users/register/). After your request has been approved, you'll see your personal API key after you login. Please keep this key to yourself!

2. `02_generate_simput.py`: Create SIMPUT files based on the previously downloaded files.

3. `03_xmm_simulation.py`: Simulate XMM-Newton for the previously created SIMPUT files. TBD: I will add at least one other satellite to choose from.

4. `04_combine_simulations.py`: **Not used right now!** I will rewrite this step to merge images from different satellites/different sensors.

### Configuring the process
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
- `num_processes`: How many processes should be run asynchrounosly. My recommendation: Your available RAM divided by 4, i. e., if you have 32GB of RAM: $32 / 4 = 8$
- `top_n`: How many cutouts to download for given simulations.
- `resolutions`: For every downloaded cutout, create images in these given resolutions

Executing the code is same for both setups:

1. Set your configuration parameters as needed (see below)
2. Choose what step you want to run (see below)
2. Run `conda run -n xmm --no-capture-output python /path/to/script -p /home/studtodorkov/simulator/config.json`

## Illustris TNG simulations `illustris_tng_image_gen.py`
For our XMM simulations we need sources to simulate (simulation input). 
In our project we are especially interested in extended sources. 
We take these extended sources from the Illustris TNG project (https://www.tng-project.org/).
This is a large cosmological hydrodynamical simulation of galaxy formation containing hundreds
of terabytes of simulated data. From this we take the most massive objects and take 
x-ray projections and x-ray slices (less realistic but contains more clearly defined structure).
Note that cutout files are relatively large (100-1000 mb) and can take a while to download, it will first download all 
the relevant cutouts before generating the images.


#### Notes on extended sources (images as simput)
To create the simput for extended sources we use fits image files. In order to have a realistic distribution we augment these images using:
- Brightness: The brightness of the source is internally defined as sigma_b. This is based on the std of 50ks background. I.e. `sigma_b = 10` will result in a brightness of 10 times the background at 50ks.
The images are used as a distribution of a given brightness. 
We determine the final brightness by taking a center cutout of the image and set this to the brightness defined by sigma_b. This behaviour can be changed in `simput/img_simputgen.py`
- Location: We augment to location by offsetting the image from the bore-axis. Since real xmm observation are usually focussed on the center of extended sources we by default offcenter the images by a small amount around the bore-sight based on a normal distribution.
- Size (zoom): We augment the size of the extended source by artificially zooming in or out.

#### Some parameters to consider:
- `home`: Change this to the data location, make sure this is the same path as used in `illustris_tng_image_gen.py` in order to use the illustris tng generated images.
- `simput_in_image_dataset`: The name of the directory containing fits images used in the image mode. By in the default workflow these will be the illustris tng images.
- `num`: The number of simputs to generate of a certain mode. For the `img` mode, if set to `-1` it will process every image
- `num_img_sample`: How many variations to generate of one image. These variations consist out of the brightness, location and zoom
- `zoom_img_range`: This the size of the object by zooming into the object.
- `sigma_b_img_range`: The brightness sample range. This is based on the std of 50ks background. 
I.e. `sigma_b = 10` will result in a brightness of 10 times the background at 50ks.
- `offset_std`: The standard deviation of the normal distribution of the offset location around the bore-sight

## XMM Simulation `xmm_simulation.py`
The XMM simulations are done using SIXTE X-ray simulation software (https://www.sternwarte.uni-erlangen.de/research/sixte/).
All the elements that make up a XMM observation are simulated separately: extended source, agn and background.
These can then in the future be combined with a detector-mask to create a realistic XMM observation.
Since this is a simulation we can also simulate observations where XMM has a higher resolution (both spatial and psf wise).

## Acknowledgements
