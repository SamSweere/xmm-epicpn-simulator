# XMM Simulation
The easiest way to run this code is to use the accompanying docker container. All the dependencies and filepaths are already pre-installed and set.


## Setup (non docker):
 - Install Heasoft (if not presently installed): https://heasarc.gsfc.nasa.gov/lheasoft/install.html
 - Install SIXTE: https://www.sternwarte.uni-erlangen.de/research/sixte/simulation.php
 - Add SIXTE to the path: `export PATH=${PATH}:/path_to_sixte_bin`
 - Download xmm instrument files: http://vospace.esac.esa.int/vospace/sh/adfa36d39939809c41fab0c6bdf3049661f6dba?dl=1 (MD5 Hash: 8018ae18bcfe0ddaa8edd353f8a5bd34)
 - Make the instruments dicrectory at: `sixte/share/sixte/instruments/`
 - Unpack the instrument files in: `sixte/share/sixte/instruments/`
 - Clone the code repository: `git clone https://github.com/SamSweere/xmm_simulation.git`
 - Open the repository: `cd xmm_simulation`
 - Create a virtual environment: `python3 -m venv xmm_simulation_venv`
 - Activate the environment: `source xmm_simulation_venv/bin/activate`
 - Install the requirements: `pip3 install -r requirements.txt`
 
## Setup Docker:
Docker is a containerised system which makes it possible create and run a whole operation system in a container.
This way all the programs are already installed. You only need to download and run the docker. 
Docker works on Linux, Windows and Mac.
- Download the docker software: 
- If not already installed, install docker: https://www.docker.com/
- Download the docker file: http://vospace.esac.esa.int/vospace/sh/a6314f6fcdeb447e4cdeb5351126df26d6aa032?dl=1 (MD5 Hash: 8018ae18bcfe0ddaa8edd353f8a5bd34)
- Load the docker file: `docker load --input xmm_sim_docker.tar`
- Create an external volume where the data is being saved: `mkdir {path_to_your_data_directory}`
- Give the folder read write access for all (docker runs as another user): `chmod a+w {path_to_your_data_directory}`
- Run the docker, to be able to acess the data we also mount an external volume: `docker run -it --mount type=bind,source={path_to_your_data_directory},target=/home/heasoft/data samsweere/sixte_xmm_heasoft:latest bash`
- Switch to the heasoft user: `su -l heasoft`
- Enable autocomplete with bash: `bash`
- Navigate to the xmm_simulation: `cd xmm_simulation`
- Run the desired files: `python {the_desired_python_file.py}`
- You can change the parameters of the files by using nano: `nano {the_desired_python_file.py}`
- Once the changes have been made in nano you can save them with `CTRL + o` followed with `enter` and exit nano with `CTRL + x`

## Data workflow
This repository contains all the code to create simulated XMM images containing extended sources, agns and background. 
All the parts can be run separately, however they might need data from previous steps. 
Every step is more elaborately explained in further sections of this readme. 
All parameters represent what I used for my research. These can be changed to test different things but might break next steps.
The workflow:
- Download and process Illustris TNG simulations: `illustris_tng_image_gen.py`
- Create simput (simulation input) files: `simput_gen.py`
- Run XMM SIXTE simulations: `xmm_simulation.py`
- (Optional) Combine simulation outputs to create an XMM like observation: `combine_simulations.py`

## Illustris TNG simulations `illustris_tng_image_gen.py`
For our XMM simulations we need sources to simulate (simulation input). 
In our project we are especially interested in extended sources. 
We take these extended sources from the Illustris TNG project (https://www.tng-project.org/).
This is a large cosmological hydrodynamical simulation of galaxy formation containing hundreds
of terabytes of simulated data. From this we take the most massive objects and take 
x-ray projections and x-ray slices (less realistic but contains more clearly defined structure).
Note that cutout files are relatively large (100-1000 mb) and can take a while to download, it will first download all 
the relevant cutouts before generating the images.

#### API-key
In order to download the Illustris TNG simulations one needs an api key. This is not included in the project.
The api key can be requested at: https://www.tng-project.org/users/register/ <br>
Once you have an api key put it in a text file at: `illustris_tng/api_key.txt`. This only has to be done once.

#### Some parameters to consider:
- `home`: Change this to the location to where you want the files to be stored
- `simulation_names`: Select which TNG simulations to use. For a quick test set this to one simulation.
- `top_n`: How many images to select from the simulations, change to a small number for testing
- `modes`: Select if you want projections, slices or both 

The code will take generate images with all the combinations of parameters. 
Thus, one TNG simulation cutout will by default generate multiple images.

## Simulation Input (simput) generation `simput_gen.py`
In order to run the XMM SIXTE simulations we need input for these simulations. 
In this project we focus on:
 - Extended sources, provided in the form of fits images, such as the Illustris TNG images
 - AGN's, generated based on real xmm observed agn distributions 
 - Background, generated based on real xmm background
 - (Optional) Test grid, a grid containing sources of the same brightness and a source on the bore-axis. Intended for use in development and testing. 

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

#### Some parameters to consider:
 - `home`: Change this to the data location, make sure this is the same path as used in `simput_gen.py` in order to use the generated simputs.
 - `instrument_dir`: SIXTE instrument location
 - `exposure`: The exposure time to simulate. The final observation will also be split into shorter observation images in steps of 10ks
 - `res_mult`: The resolution (both spatial and psf) multiplier. `1` matches real xmm observations, `2` and `4` increase the spatial resolution by 2 and 4 and decrease the psf by 2 and 4 respectively.
 - `amount`: The amount of simulations to do per mode. `-1` will run the simulation on all the simputs of the mode

## Combine Simulations `combine_observations.py`
The simulator will produce all separate parts of the xmm observation. To create a realistic xmm observation we need to combine these simulated parts.
We have the option to select the separate parts, finally the image is multiplied with the detector mask

#### Some parameters to consider:
- `home`: Change this to the data location, make sure this is the same path as used in `xmm_simulation.py` in order to use the simulated observation parts.
- `num`: The number of images to combine from this mode, `-1` will combine everything image in this mode
- `agn`: Option to add AGNS
- `background`: Option to add background
- `sample_n`: The number of agns and backgrounds to sample for one observation
- `res_mult`: The resolution multiples to use. It will find the same simput simulation files for the multiple resolutions
- `exposure`: The exposures to combine
