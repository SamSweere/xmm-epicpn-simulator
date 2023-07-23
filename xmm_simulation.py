import os
import shutil
from pathlib import Path

from sixte.simulator import run_simulation_mode
from utils import log

config = {
    'docker': False,
    'home': os.path.join(os.path.expanduser("~"), 'Documents/ESA/data/sim'),
    # The home directory when not running inside docker
    'instrument_dir': os.path.join(os.path.expanduser("~"),
                                   'Documents/ESA/programs/sixte/share/sixte/instruments/xmm/epicpn'),
    # The location of the SIXTE instruments directory, no need to change when running in docker
    'debug': False,
    'verbose': False,
    'gb_per_process': 0.5,
    'num_processes': 0,  # 0 will use all threads of the cpu or limits it by the memory capability
    'exposure': [100000],
    'data_dir': None,
    'dataset_dir': None,
    'run_dir': None,
    'simput_base_path': None,
    'res_mult': [1, 2],  # Options: [1,2,4]
    'mode': ['img', 'agn', 'background', 'test_grid'],  # Options: 'img', 'agn', 'background', 'test_grid', 'exposure'
    'amount': [-1, -1, 10000, -1],  # -1 will generate all the
    'dataset_name': "xmm_sim_dataset",
    'pi0': 500,
    'pi1': 2000,
    'simput_order': 'normal',
    # normal (fron to back), reversed (back to front), random.  Order of processing the simputs,
    # random helps with collisions
    'simulate_separate_ccds': False,
    # Option to instead of simulating one big ccd, to simulate separate ccd's
    # #TODO: the 2x and 4x allignment of the ccd's is incorrect!
}

# Set the data and sas dir depending on docker or local use
if config['docker']:
    config['data_dir'] = '/home/heasoft/data/sim'
    config['run_dir'] = '/home/heasoft/tmp'
    config['instrument_dir'] = '/home/heasoft/sixte/share/sixte/instruments/xmm/epicpn'
else:
    # Running on laptop/pc
    config['data_dir'] = config['home']
    config['run_dir'] = os.path.join(config['data_dir'], 'tmp')

config['simput_base_path'] = os.path.join(config['data_dir'], 'simput')
config['dataset_dir'] = os.path.join(config['data_dir'], config['dataset_name'])

if __name__ == '__main__':
    # Setup logging
    log.setup_logging(config, prefix="sixte_sim")

    # Create dataset dir if it does not yet exist
    if not os.path.exists(config['dataset_dir']):
        os.makedirs(config['dataset_dir'])

    # Copy the masks to the dataset dir
    mask_dir = os.path.join(config['dataset_dir'], 'detector_mask')
    if not os.path.exists(mask_dir):
        os.makedirs(mask_dir)

    file_path = os.path.dirname(__file__)
    for res_mult in config['res_mult']:
        res_maks_dir = os.path.join(mask_dir, f'{res_mult}x')
        if not os.path.exists(res_maks_dir):
            os.makedirs(res_maks_dir)

        mask_file_name = f'pn_mask_500_2000_detxy_{res_mult}x.ds'
        src_mask_path = os.path.join(file_path, 'detector_masks', mask_file_name)
        dst_mask_path = os.path.join(res_maks_dir, mask_file_name)
        if not os.path.exists(dst_mask_path):
            shutil.copyfile(src_mask_path, dst_mask_path)

    mode_list = config['mode']
    amount_list = config['amount']

    if not (len(mode_list) == len(amount_list)):
        m = f'the number of modes ({len(mode_list)}) should be equal to the number of amounts ({len(amount_list)}) '
        log.elog(m)
        raise AssertionError(m)

    # Create tmp run_dir
    run_dir = Path(config['run_dir'])
    run_dir.mkdir(parents=True, exist_ok=True)

    run_simulation_mode(mode_list, amount_list, run_dir, config=config)

    for mode, amount in zip(mode_list, amount_list):
        # Process the mode with the settings
        print(f"Running {amount} simulations with mode: {mode}")
        run_simulation_mode(mode=mode, amount=amount, run_dir=run_dir, config=config)

    if not config['debug']:
        # Remove the tmp run_dir
        shutil.rmtree(run_dir)
