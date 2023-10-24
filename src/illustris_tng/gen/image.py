import json
from multiprocessing.pool import Pool
from pathlib import Path
from typing import List, Dict, Union, Any

from loguru import logger

from src.illustris_tng.data import get, get_cutouts
from src.illustris_tng.gen.fits import cutout_to_xray_fits
from src.xmm_utils.multiprocessing import get_num_processes


def run(
        path_to_cfg: Path,
        api_key: str,
        cloudy_emissivity_root: Path
):
    with open(path_to_cfg, "r") as file:
        cfg: Dict[str, dict] = json.load(file)
    env_cfg: Dict[str, Any] = cfg["environment"]
    mp_cfg: Dict[str, Any] = cfg["multiprocessing"]
    illustris_cfg: Dict[str, Any] = cfg["illustris"]

    debug = env_cfg["debug"]

    log_dir = Path(env_cfg["log_directory"])
    log_dir.mkdir(parents=True, exist_ok=True)
    log_level = "DEBUG" if env_cfg["debug"] else "INFO"
    logger.add(f"{(log_dir / '01_download_files_{time}.log')}", enqueue=True,
               format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}", level=log_level)

    if not debug:
        logger.info(f"Since 'debug' is set to 'false' the download will be run asynchronously.")

    working_dir = Path(env_cfg["working_directory"]).expanduser()
    cutout_dir = working_dir / "cutout_data"
    dataset_dir = working_dir / illustris_cfg["dataset_dir"]
    dataset_dir.mkdir(parents=True, exist_ok=True)

    base_url = "http://www.tng-project.org/api/"
    headers = {"api-key": api_key}

    logger.info("Getting snaps...")
    simulations: List[Dict[str, Union[str, int]]] = get(base_url, headers=headers, cutout_datafolder=cutout_dir)[
        "simulations"]
    simulations_to_use: List[Dict[str, Union[str, int]]] = []
    for simulation in simulations:
        if simulation["name"] in illustris_cfg["simulation_names"]:
            simulations_to_use.append(simulation)

    subs_r = []  # Save the sub results in one list

    for simulation in simulations_to_use:
        name = simulation["name"]
        snapshots_url: str = get(path=simulation["url"], headers=headers, cutout_datafolder=cutout_dir)["snapshots"]
        snaps: List[Dict[str, Any]] = get(path=snapshots_url, headers=headers, cutout_datafolder=cutout_dir)
        snaps = [snaps[z] for z in illustris_cfg['snapshot_nums']]

        subs_sim_r = []
        logger.info(f"\tGetting subs for: {name}")
        # There are 100 snapshots, the last one corresponds to z=0
        for z, snap in zip(illustris_cfg['snapshot_nums'], snaps):
            logger.info(f"\tsnapshot num (z): {z}")
            subhalos = get(snap['url'], headers=headers, cutout_datafolder=cutout_dir)["subhalos"]

            # request and inspect most massive 100 subhalos that are central (primary_flag = 1)
            # primary_flag = 1 indicates that this is the central (i.e. most massive, or "primary") subhalo of this FoF halo.

            subs_sim_z_r = get(subhalos, headers=headers, cutout_datafolder=cutout_dir,
                               params={'limit': illustris_cfg['top_n'], 'primary_flag': 1,
                                       'order_by': '-mass_gas'})["results"]  # 'order_by':'-mass_dm'
            subs_sim_r.extend(get_cutouts(subs=subs_sim_z_r, headers=headers, cutout_datafolder=cutout_dir,
                                          fail_on_error=env_cfg["fail_on_error"]))

        subs_r += subs_sim_r
        logger.info(f"\tNumber of subs loaded: {len(subs_sim_r)}")

    logger.info(f"Total subs loaded: {len(subs_r)}")

    for i, mode in enumerate(illustris_cfg["modes"].keys()):
        if mode not in ("proj", "slice"):
            raise ValueError(f"Expected either 'proj' or 'slice' at index {i} but got mode '{mode}'!")

    with Pool(get_num_processes(mp_conf=mp_cfg)) as pool:
        for res in subs_r:
            sub = res["sub"]
            cutout = res["cutout"]
            # Generate the fits file
            # We can change the normal to rotate the image
            cutout_args = (cutout, dataset_dir, sub, illustris_cfg['modes'], cloudy_emissivity_root,
                           illustris_cfg['emin'], illustris_cfg['emax'], illustris_cfg['width'],
                           illustris_cfg['resolutions'], illustris_cfg['redshift'], illustris_cfg['overwrite'])
            if debug:
                pool.apply(cutout_to_xray_fits, cutout_args)
            else:
                pool.apply_async(cutout_to_xray_fits, cutout_args)
        pool.close()
        pool.join()
    logger.info("Done!")
