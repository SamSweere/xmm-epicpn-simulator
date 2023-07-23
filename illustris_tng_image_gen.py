import os
from pathlib import Path
from typing import List, Dict, Union, Any

from illustris_tng.fits_gen import cutout_to_xray_fits
from src.illustris_tng.download_data import get, get_and_save_cutouts
from utils.multiprocessing import run_apply_async_multiprocessing

config = {
    'docker': False,
    'home': os.path.join(os.path.expanduser("~"), 'Documents/ESA/data/sim'),
    # Change this if not running in docker mode
    'gb_per_process': 2.0,
    'num_processes': 0,  # 0 is as many as possible within the cpu/memory capabilities
    'dataset_name': "tng300_1_2_3_z99_2048",  # The name of directory the files will be saved in
    'simulation_names': ['TNG300-3'],
    # 'TNG300-2', 'TNG300-3'],  # Options: ['TNG50-1', 'TNG50-2', 'TNG50-3', 'TNG50-4',
    # 'TNG100-1', 'TNG100-2', 'TNG100-3', 'TNG300-1', 'TNG300-2', 'TNG300-3'],
    # Illlustris TNG300 is the largest simulation (300Mpc box) and has 3 resolutions and 3 sub boxes of full physics.
    'snapshot_nums': [99],  # The tng simulations have snapshots at different times, 99 is the last snapshot,
    # one can add multiple time snapshots with a comma. I.e. [1, 50, 99]
    'top_n': 100,  # 100 #Top number of each simulation
    # The cutout_to_simput is quite time_consuming therefore we multi-process it
    'emin': 0.5,
    'emax': 2.0,
    'redshift': 0.05,  # The redshift to the object.
    'overwrite': False,
    'create_preview': False,  # Create preview images
    'width': [1000, 4000],  # How much to zoom by changing the width of the image
    'resolutions': [2048],
    'modes': ['proj', 'slice'],  # ['proj', 'slice']
    # Take the normals from all directions
    'proj_normals': [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
    # Normals do not work for slices, thus make them axes
    'slice_axes': ['x', 'y', 'z'],
}

# Change the run directory to illustris_tng such that the illustris code finde the cloudy_emissivity_v2.h5 file
# root_dir = os.path.dirname(__file__)
# run_dir = tng_api_key_path = os.path.join(root_dir, 'illustris_tng')
# os.chdir(run_dir)

root_dir = Path.cwd()

# Load the illustris tng api_key
# TODO Make this customizable
tng_api_key_path = root_dir / "api_key.txt"

if tng_api_key_path.exists():
    # TODO this is ugly
    with open(tng_api_key_path, "r") as txt:
        api_key = txt.readline().replace("\n", "")
else:
    raise FileNotFoundError

if config['docker']:
    dataset_base_dir: Path = '/home/heasoft/data/sim'
else:
    dataset_base_dir: Path = Path.home() / "testing"

cutout_dir = dataset_base_dir / "cutout_data"
dataset_dir = dataset_base_dir / config["dataset_name"]
dataset_dir.mkdir(parents=True, exist_ok=True)

baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key": api_key}

print("Getting snaps...")

simulations: List[Dict[str, Union[str, int]]] = get(baseUrl, headers=headers, cutout_datafolder=cutout_dir)[
    "simulations"]
simulations_to_use: List[Dict[str, Union[str, int]]] = []
for simulation in simulations:
    if simulation["name"] in config["simulation_names"]:
        simulations_to_use.append(simulation)

subs_r = []  # Save the sub results in one list
for simulation in simulations_to_use:
    name = simulation["name"]
    simulation: Dict[str, Any] = get(path=simulation["url"], headers=headers, cutout_datafolder=cutout_dir)
    snaps: List[Dict[str, Any]] = get(path=simulation["snapshots"], headers=headers, cutout_datafolder=cutout_dir)

    subs_sim_r = []
    print(f"Getting subs of {name}")
    # There are 100 snapshots, the last one corresponds to z=0
    for z in config['snapshot_nums']:
        print(f"snapshot num (z): {z}")
        snap = get(snaps[z]['url'], headers=headers, cutout_datafolder=cutout_dir)

        # request and inspect most massive 100 subhalos that are central (primary_flag = 1)
        # primary_flag = 1 indicates that this is the central (i.e. most massive, or "primary") subhalo of this FoF halo.

        subs_sim_z_r = get(snap['subhalos'], headers=headers, cutout_datafolder=cutout_dir,
                           params={'limit': config['top_n'], 'primary_flag': 1,
                                   'order_by': '-mass_gas'})  # 'order_by':'-mass_dm'
        subs_sim_r += subs_sim_z_r['results']

    subs_r += subs_sim_r
    print(f"Number of subs loaded for {name}", len(subs_sim_r))

print(f"Total subs loaded", len(subs_r))

# Download all the cutouts
sc = get_and_save_cutouts(subs=subs_r, headers=headers, cutout_datafolder=cutout_dir)

argument_list = []

for sub, cutout in sc:
    # Generate the fits file
    # We can change the normal to rotate the image
    for mode in config['modes']:
        if mode == "proj":
            normals = config['proj_normals']
        elif mode == "slice":
            normals = config['slice_axes']
        else:
            raise ValueError(f"mode {mode} not suported")

        for width in config['width']:
            for normal in normals:
                for resolution in config['resolutions']:
                    args = (cutout, cutout_dir, dataset_dir, sub, mode, config['emin'], config['emax'],
                            normal, width, resolution, config['redshift'], config['overwrite'],
                            config['create_preview'])
                    argument_list.append(args)

        # For debugging, this does not use multiprocessing
        # cutout_to_xray_fits(cutout, cutout_datafolder, dataset_dir, sub, mode, emin, emax,
        #                                             normal, width, resolution, redshift, overwrite,
        #                                             create_preview)

run_apply_async_multiprocessing(cutout_to_xray_fits, argument_list, num_processes=config['num_processes'],
                                gb_per_process=config['gb_per_process'])
