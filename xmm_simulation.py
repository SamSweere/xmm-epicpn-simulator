import json
from pathlib import Path
from typing import Dict, Union

from src.sixte.simulator import run_simulation_modes


# config = {
#     'data_dir': None,
#     'dataset_dir': None,
#     'run_dir': None,
#     'simput_base_path': None,
#     'dataset_name': "xmm_sim_dataset",
#     # 'pi0': 500,
#     # 'pi1': 2000,
#     'simput_order': 'normal',
#     # normal (fron to back), reversed (back to front), random.  Order of processing the simputs,
#     # random helps with collisions
#     # Option to instead of simulating one big ccd.py, to simulate separate ccd.py's
#     # #TODO: the 2x and 4x allignment of the ccd.py's is incorrect!
# }
#
# # Set the data and sas dir depending on docker or local use
# if config['docker']:
#     config['data_dir'] = '/home/heasoft/data/sim'
#     config['run_dir'] = '/home/heasoft/tmp'
#     config['instrument_dir'] = '/home/heasoft/sixte/share/sixte/instruments/xmm/epicpn'
# else:
#     # Running on laptop/pc
#     config['data_dir'] = config['home']
#     config['run_dir'] = os.path.join(config['data_dir'], 'tmp')
#
# config['simput_base_path'] = os.path.join(config['data_dir'], 'simput')
# config['dataset_dir'] = os.path.join(config['data_dir'], config['dataset_name'])

# if __name__ == '__main__':
#     # Setup logging
#
#     # Create dataset dir if it does not yet exist
#     if not os.path.exists(config['dataset_dir']):
#         os.makedirs(config['dataset_dir'])
#
#     # Copy the masks to the dataset dir
#     mask_dir = os.path.join(config['dataset_dir'], 'detector_mask')
#     if not os.path.exists(mask_dir):
#         os.makedirs(mask_dir)
#
#     file_path = os.path.dirname(__file__)
#     for res_mult in config['res_mult']:
#         res_maks_dir = os.path.join(mask_dir, f'{res_mult}x')
#         if not os.path.exists(res_maks_dir):
#             os.makedirs(res_maks_dir)
#
#         mask_file_name = f'pn_mask_500_2000_detxy_{res_mult}x.ds'
#         src_mask_path = os.path.join(file_path, 'res/detector_masks', mask_file_name)
#         dst_mask_path = os.path.join(res_maks_dir, mask_file_name)
#         if not os.path.exists(dst_mask_path):
#             shutil.copyfile(src_mask_path, dst_mask_path)
#
#     mode_list = config['mode']
#     amount_list = config['amount']
#
#     # Create tmp run_dir
#     run_dir = Path(config['run_dir'])
#     run_dir.mkdir(parents=True, exist_ok=True)
#
#     for mode, amount in zip(mode_list, amount_list):
#         # Process the mode with the settings
#         print(f"Running {amount} simulations with mode: {mode}")
#         run_simulation_modes(mode=mode, amount=amount, run_dir=run_dir, config=config)
#
#     if not config['debug']:
#         # Remove the tmp run_dir
#         shutil.rmtree(run_dir)


def run(cfg: Union[Path, Dict[str, dict]]) -> None:
    if isinstance(cfg, Path):
        with open(cfg, "r") as f:
            cfg: Dict[str, dict] = json.load(f)

    env_cfg = cfg["environment"]
    simput_cfg = cfg["simput"]
    xmm_cfg = cfg["xmm"]

    working_directory = Path(env_cfg["working_directory"]).expanduser()

    dataset_dir = working_directory / "xmm_sim_dataset"
    dataset_dir.mkdir(parents=True, exist_ok=True)

    simput_path = working_directory / "simput"
    if not simput_path.exists():
        raise NotADirectoryError(f"Simput directory '{simput_path.resolve()}' not found! "
                                 f"Please make sure that the simput files have been generated and that "
                                 f"the working directory is set correctly.")

    mode_amount_dict: Dict[str, int] = xmm_cfg["mode"]
    run_simulation_modes(mode_amount_dict, working_directory, simput_path, simput_cfg["order"], cfg)


if __name__ == '__main__':
    run(Path("/home/bojantodorkov/Projects/xmm-epicpn-simulator/cfg/xmm.json"))
