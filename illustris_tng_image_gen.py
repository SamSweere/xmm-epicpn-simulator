import os

from illustris_tng.download_data import get, get_cutout_from_sub, get_and_save_cutouts
from illustris_tng.fits_gen import cutout_to_xray_fits
from utils.multiprocessing import run_apply_async_multiprocessing

config = {
    'docker': True,
    'home': os.path.join(os.path.expanduser("~"), 'Documents/ESA/data/sim'),  # Change this if not running in docker mode
    'gb_per_process': 2.0,
    'num_processes': 0,  # 0 is as many as possible within the cpu/memory capabilities
    'dataset_name': "tng300_1_2_3_z99_2048",  # The name of directory the files will be saved in
    'simulation_names': ['TNG300-1', 'TNG300-2', 'TNG300-3'],  # Options: ['TNG50-1', 'TNG50-2', 'TNG50-3', 'TNG50-4',
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

if config['docker']:
    dataset_base_dir = '/home/heasoft/data/sim'
else:
    dataset_base_dir = config['home']

cutout_dir = os.path.join(dataset_base_dir, "cutout_data")
working_dir = os.path.join(dataset_base_dir, config['dataset_name'])

# Change the run directory to illustris_tng such that the illustris code finde the cloudy_emissivity_v2.h5 file
root_dir = os.path.dirname(__file__)
run_dir = tng_api_key_path = os.path.join(root_dir, 'illustris_tng')
os.chdir(run_dir)

# Load the illustris tng api_key
tng_api_key_path = os.path.join(root_dir, 'api_key.txt')

if not os.path.exists(tng_api_key_path):
    raise FileNotFoundError(f"No Illustris Tng api key found at: {tng_api_key_path}. "
                            f"Please add the illustris tng api key to a text file at: 'illustris_tng/api_key.txt'")
else:
    infile = open(tng_api_key_path, 'r')
    tng_api_key = infile.readline().replace('\n', '')

dataset_dir = os.path.join(dataset_base_dir, config['dataset_name'])

if not os.path.exists(dataset_dir):
    # If the folder does not exist make it
    os.makedirs(dataset_dir)
    print("Created dataset_dir:", dataset_dir)

cutout_datafolder = cutout_dir

baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key": tng_api_key}

print("Getting snaps...")

r = get(baseUrl, headers=headers, cutout_datafolder=cutout_datafolder)
names = [sim['name'] for sim in r['simulations']]

subs_r = []  # Save the sub results in one list
for name in config['simulation_names']:
    i = names.index(name)
    sim = get(r['simulations'][i]['url'], headers=headers, cutout_datafolder=cutout_datafolder)
    snaps = get(sim['snapshots'], headers=headers, cutout_datafolder=cutout_datafolder)

    subs_sim_r = []
    print("Getting subs of", name)
    # There are 100 snapshots, the last one corresponds to z=0
    for z in config['snapshot_nums']:
        print("snapshot num (z):", z)
        snap = get(snaps[z]['url'], headers=headers, cutout_datafolder=cutout_datafolder)

        # request and inspect most massive 100 subhalos that are central (primary_flag = 1)
        # primary_flag = 1 indicates that this is the central (i.e. most massive, or "primary") subhalo of this FoF halo.

        subs_sim_z_r = get(snap['subhalos'], headers=headers, cutout_datafolder=cutout_datafolder,
                           params={'limit': config['top_n'], 'primary_flag': 1,
                                   'order_by': '-mass_gas'})  # 'order_by':'-mass_dm'
        subs_sim_r += subs_sim_z_r['results']

    subs_r += subs_sim_r
    print(f"Number of subs loaded for {name}", len(subs_sim_r))

print(f"Total subs loaded", len(subs_r))

# Download all the cutouts
get_and_save_cutouts(subs=subs_r, headers=headers, cutout_datafolder=cutout_datafolder)

argument_list = []

for sub_num in range(len(subs_r)):
    sub = get(subs_r[sub_num]['url'], headers=headers, cutout_datafolder=cutout_datafolder)
    cutout = get_cutout_from_sub(sub=sub,
                                 headers=headers,
                                 cutout_datafolder=cutout_datafolder)  # If the cutout is aldready downloaded it loads it, otherwise it will download and save the cutout

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
                    args = (cutout, cutout_datafolder, dataset_dir, sub, mode, config['emin'], config['emax'],
                            normal, width, resolution, config['redshift'], config['overwrite'],
                            config['create_preview'])
                    argument_list.append(args)

                    # For debugging, this does not use multiprocessing
                    # cutout_to_xray_fits(cutout, cutout_datafolder, dataset_dir, sub, mode, emin, emax,
                    #                                             normal, width, resolution, redshift, overwrite,
                    #                                             create_preview)

run_apply_async_multiprocessing(cutout_to_xray_fits, argument_list, num_processes=config['num_processes'],
                                gb_per_process=config['gb_per_process'])
