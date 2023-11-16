import json
from argparse import ArgumentParser
from datetime import datetime, timedelta
from functools import partial
from pathlib import Path
from typing import Any, Dict, List

from loguru import logger

from src.illustris_tng.download_data import get_available_simulations, get_cutouts, get_subhalos
from src.illustris_tng.fits import cutout_to_xray_fits
from src.xmm_utils.multiprocessing import get_pool
from src.xmm_utils.run_utils import configure_logger, create_dirs, handle_error

logger.remove()


def run(
        path_to_cfg: Path,
        api_key: str,
        cloudy_emissivity_root: Path
):
    starttime = datetime.now()
    with open(path_to_cfg, "r") as file:
        cfg: Dict[str, dict] = json.load(file)
    env_cfg: Dict[str, Any] = cfg["environment"]
    mp_cfg: Dict[str, Any] = cfg["multiprocessing"]
    illustris_cfg: Dict[str, Any] = cfg["illustris"]

    debug = env_cfg["debug"]

    for i, mode in enumerate(illustris_cfg["modes"].keys()):
        if mode not in ("proj", "slice"):
            raise ValueError(f"Expected either 'proj' or 'slice' at index {i} but got mode '{mode}'!")

    # Create all directories
    log_dir = Path(env_cfg["log_directory"]).expanduser()
    working_dir = Path(env_cfg["working_directory"]).expanduser()
    cutout_dir = working_dir / "cutout_data"
    dataset_dir = working_dir / illustris_cfg["dataset_dir"]
    create_dirs([log_dir, working_dir, cutout_dir, dataset_dir])

    configure_logger(log_dir=log_dir, log_name="01_download_files.log", enqueue=True, debug=debug,
                     rotation=timedelta(hours=1))

    if not debug:
        logger.info(f"Since 'debug' is set to 'false' the download will be run asynchronously.")

    simulations: List[dict] = get_available_simulations(api_key=api_key)
    filtered_simulations: List[str] = []
    for simulation in simulations:
        if simulation["name"] in illustris_cfg["simulation_names"]:
            filtered_simulations.append(simulation["url"])
    logger.info(f"Found {len(filtered_simulations)} available simulations. "
                f"({[sim for sim in filtered_simulations]})")

    subhalos = []
    for simulation in filtered_simulations:
        for snapshot_num in illustris_cfg['snapshot_nums']:
            # request and inspect most massive 100 subhalos that are central (primary_flag = 1)
            # primary_flag = 1 indicates that this is the central (i.e. most massive, or "primary") subhalo of
            # this FoF halo.
            subhalos.extend(get_subhalos(api_key=api_key,
                                         simulation_url=simulation,
                                         snapshot_num=snapshot_num,
                                         params={'limit': illustris_cfg['top_n'], 'primary_flag': 1,
                                                 'order_by': '-mass_gas'}))

    with get_pool(mp_conf=mp_cfg) as pool:
        mp_apply = pool.apply if debug else partial(pool.apply_async, error_callback=handle_error)
        logger.info("START\tDownloading/getting already downloaded cutouts.")
        for subhalo in subhalos:
            arguments = (subhalo, api_key, cutout_dir, env_cfg["fail_on_error"])
            mp_apply(get_cutouts, arguments)

        pool.close()
        pool.join()
        logger.info("DONE\tDownloading/getting cutouts.")

    with get_pool(mp_conf=mp_cfg) as pool:
        mp_apply = pool.apply if debug else partial(pool.apply_async, error_callback=handle_error)
        logger.info("START\tGenerating FITS from cutouts.")
        cutouts = cutout_dir.glob("*.hdf5")
        for cutout_path in cutouts:
            pos_x, pos_yz = cutout_path.stem.split("_x_")[1].split("_y_")
            pos_y, pos_z = pos_yz.split("_z_")
            sub = {
                "pos_x": float(pos_x),
                "pos_y": float(pos_y),
                "pos_z": float(pos_z)
            }

            arguments = (cutout_path, dataset_dir, sub, illustris_cfg['modes'], cloudy_emissivity_root,
                         illustris_cfg['emin'], illustris_cfg['emax'], illustris_cfg['width'],
                         illustris_cfg['resolutions'], illustris_cfg['redshift'], illustris_cfg['overwrite'])
            mp_apply(cutout_to_xray_fits, arguments)
        pool.close()
        pool.join()
        logger.info("DONE\tGeneration of FITS done.")
    endtime = datetime.now()
    logger.info(f"Duration: {endtime - starttime}")


if __name__ == '__main__':
    parser = ArgumentParser(prog="", description="")
    parser.add_argument("-k", "--api_key", type=str, required=True,
                        help="IllustrisTNG API key. If you don't have one, create an account at "
                             "https://www.tng-project.org/data/")
    parser.add_argument("-p", "--config_path", type=Path, required=True, help="Path to config file.")
    parser.add_argument("-c", "--cloudy_emissivity", type=Path, required=True,
                        help="Path to root directory of cloudy_emissivity_v2.h5")

    args = parser.parse_args()
    run(path_to_cfg=args.config_path, api_key=args.api_key, cloudy_emissivity_root=args.cloudy_emissivity)
