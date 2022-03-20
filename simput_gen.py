import os
import shutil

import numpy as np
from tqdm import tqdm

from simput.generation import generate_simput
from utils import log
from utils.multiprocessing import get_multiprocessing_pool

config = {
    'docker': True,
    'debug': False,
    'verbose': False,
    'home': os.path.join(os.path.expanduser("~"), 'Documents/ESA/data/sim'),  # Change this if not running in docker mode
    'gb_per_process': 0.1,
    'num_processes': 0,
    'data_dir': None,
    'run_dir': None,
    'simput_config': None,
}

simput_config = {
    'mode': ['img', 'agn', 'background', 'test_grid', 'exposure_map'],
    'num': [-1, 20000, 1, -1, 0],  # -1 for all
    'zoom_img_range': [1, 2],  # Sample range
    'sigma_b_img_range': [5, 50],  # Sample range
    'offset_std': 0.05,  # Offset standard deviation in percentage of the image
    'num_img_sample': 5,
    'simput_in_image_dataset': "tng300_1_2_3_z99_2048",  # The name of the image folder to be used in the img mode,
    # the images need to be in fits format
    'sim_in_root_path': None,
    'sim_in_dataset_path': None,
    'simput_dir': None,
}

# Set the data and sas dir depending on docker or local use
if config['docker']:
    config['data_dir'] = '/home/heasoft/data/sim'

    config['run_dir'] = '/home/heasoft/tmp'
    simput_config['sim_in_root_path'] = config['data_dir']
    simput_config['simput_dir'] = os.path.join(config['data_dir'], 'simput')
    simput_config['simput_base_path'] = os.path.join(config['data_dir'], 'simput')
    simput_config['sim_in_root_path'] = config['data_dir']
    simput_config['sim_in_dataset_path'] = os.path.join(simput_config['sim_in_root_path'],
                                                        simput_config['simput_in_image_dataset'])


else:
    # Running on laptop/pc
    config['data_dir'] = config['home']
    config['run_dir'] = os.path.join(config['data_dir'], 'tmp')
    simput_config['simput_base_path'] = os.path.join(config['data_dir'], 'simput')
    simput_config['sim_in_root_path'] = config['data_dir']
    simput_config['sim_in_dataset_path'] = os.path.join(simput_config['sim_in_root_path'],
                                                        simput_config['simput_in_image_dataset'])
    simput_config['simput_dir'] = os.path.join(config['data_dir'], 'simput')

config['simput_config'] = simput_config

if __name__ == '__main__':
    # Setup logging
    log.setup_logging(config, prefix="simput_gen")

    # Create the simput dir if it does not exits
    if not os.path.exists(config['simput_config']['simput_dir']):
        os.makedirs(config['simput_config']['simput_dir'])

    # Create tmp run_dir
    run_dir = config['run_dir']
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)

    # Add all modes to the modes queue to be run
    modes = []
    for i in range(len(config['simput_config']['mode'])):
        mode = config['simput_config']['mode'][i]
        num = int(config['simput_config']['num'][i])
        if num == 0:
            continue

        img_settings = None

        # Determine the filepath for this mode, create the directory if it does not exist
        mode_file_path = os.path.join(simput_config['simput_dir'], mode)
        if not os.path.exists(mode_file_path):
            os.makedirs(mode_file_path)

        # List the already generated files
        generated_files = os.listdir(mode_file_path)

        if mode == 'img':
            # Generate the img simputs on the fly
            infiles = []
            for file in os.listdir(config['simput_config']['sim_in_dataset_path']):
                if file.endswith(".fits"):
                    infiles.append(os.path.join(config['simput_config']['sim_in_dataset_path'], file))

            if num > 0 and num < len(infiles):
                infiles = infiles[:num]

            for infile in tqdm(infiles):
                sample_num = config['simput_config']['num_img_sample']
                # Check how many of this image is already generated
                base_infilename = os.path.basename(infile).replace(".fits", "")

                for gen_file in generated_files:
                    if base_infilename in gen_file:
                        # Base_infilename in gen file, sample less
                        sample_num -= 1
                        print(
                            f"Found generated file with same basename. Reducing the number of samples for this file to {sample_num}")

                        if sample_num == 0:
                            break
                if sample_num == 0:
                    continue

                zoom_range = config['simput_config']['zoom_img_range']
                zoom_list = np.round(np.random.uniform(low=zoom_range[0], high=zoom_range[1], size=sample_num), 2)
                sigma_b_range = config['simput_config']['sigma_b_img_range']
                sigma_b_list = np.round(np.random.uniform(low=sigma_b_range[0], high=sigma_b_range[1], size=sample_num),
                                        2)

                for zoom, sigma_b in zip(zoom_list, sigma_b_list):
                    offset_std = simput_config['offset_std']
                    offset_x = round(np.random.normal(-offset_std, offset_std), 2)
                    offset_y = round(np.random.normal(-offset_std, offset_std), 2)

                    img_settings = {
                        'img_path': infile,
                        'zoom': zoom,
                        'sigma_b': sigma_b,
                        'offset_x': offset_x,
                        'offset_y': offset_y,
                    }

                    modes.append((mode, img_settings))
        else:
            for n in range(num):
                modes.append((mode, img_settings))

    # The args needed to generate a simput file
    simput_dir = config['simput_config']['simput_dir']
    debug = config['debug']
    verbose = config['verbose']

    # Run the processes parallel with a progress bar
    pool = get_multiprocessing_pool(num_processes=config['num_processes'], gb_per_process=config['gb_per_process'])

    jobs = []
    for mode, img_settings in modes:
        # For debugging:
        # generate_simput(mode, img_settings, run_dir, simput_config, debug, verbose)
        #
        # Create the output folders to not get conflicts
        # Create the simput dir for the specific mode
        output_file_dir = os.path.join(simput_dir, mode)
        if not os.path.exists(output_file_dir):
            os.makedirs(output_file_dir)

        job = pool.apply_async(generate_simput,
                               args=(mode, img_settings, run_dir, output_file_dir, simput_config, debug, verbose))
        jobs.append(job)

    pool.close()
    result_list_tqdm = []
    for job in tqdm(jobs):
        result_list_tqdm.append(job.get())

    if not config['debug']:
        # Remove the tmp run_dir
        shutil.rmtree(config['run_dir'])
