# XMM EPIC-pn Simulator
This repository hosts simulation code for generating XMM-Newton EPIC-pn observations, utilizing the [SIXTE](https://www.sternwarte.uni-erlangen.de/research/sixte/) and [SIMPUT](http://www.sternwarte.uni-erlangen.de/git.public/simput.git/) frameworks.

The simulations are tailored to produce training datasets for deep learning algorithms aimed at super-resolution enhancement and noise reduction of XMM-Newton EPIC-pn data. For details on the deep learning implementation, refer to the [xmm-superres-denoise](https://github.com/SamSweere/xmm-superres-denoise) repository. The findings from this research are documented in the paper _Deep Learning-Based Super-Resolution and De-Noising for XMM-Newton Images_, published in [MNRAS, 517, 4054 (2022)](https://doi.org/10.1093/mnras/stac2437). Please cite this publication when using the simulation code for your research.

## Installation
Given the complex dependencies and the need for external software, the code is best run in a Docker container. See the [docs/installation.md](docs/installation.md) for detailed instructions on installing the necessary software, optionally using Docker.

## Configuring the code
The good thing: You'll need to fill out `config.toml` only once! Every step relies on this configuration file and everything will be done accordingly. This file is divided into `environment`, `energy`, `download`, `simput` and `simulation`:

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
