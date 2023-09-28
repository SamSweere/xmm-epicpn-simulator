import json
import pickle
from pathlib import Path
from typing import List, Dict, Union, Any

from src.illustris_tng.data import get, get_cutouts
from src.illustris_tng.data.data_handling import save_cutout
from src.illustris_tng.gen.fits import cutout_to_xray_fits
from src.xmm_utils.multiprocessing import mp_run


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

    working_dir = Path(env_cfg["working_directory"]).expanduser()
    cutout_dir = working_dir / "cutout_data"
    dataset_dir = working_dir / illustris_cfg["dataset_dir"]
    dataset_dir.mkdir(parents=True, exist_ok=True)

    base_url = "http://www.tng-project.org/api/"
    headers = {"api-key": api_key}

    print("Getting snaps...")
    simulations: List[Dict[str, Union[str, int]]] = get(base_url, headers=headers, cutout_datafolder=cutout_dir)[
        "simulations"]
    simulations_to_use: List[Dict[str, Union[str, int]]] = []
    for simulation in simulations:
        if simulation["name"] in illustris_cfg["simulation_names"]:
            simulations_to_use.append(simulation)

    subs_r = []  # Save the sub results in one list

    for simulation in simulations_to_use:
        name = simulation["name"]
        simulation: Dict[str, Any] = get(path=simulation["url"], headers=headers, cutout_datafolder=cutout_dir)
        snaps: List[Dict[str, Any]] = get(path=simulation["snapshots"], headers=headers, cutout_datafolder=cutout_dir)

        subs_sim_r = []
        print(f"Getting subs of {name}")
        # There are 100 snapshots, the last one corresponds to z=0
        for z in illustris_cfg['snapshot_nums']:
            print(f"snapshot num (z): {z}")
            snap = get(snaps[z]['url'], headers=headers, cutout_datafolder=cutout_dir)

            # request and inspect most massive 100 subhalos that are central (primary_flag = 1)
            # primary_flag = 1 indicates that this is the central (i.e. most massive, or "primary") subhalo of this FoF halo.

            subs_sim_z_r = get(snap['subhalos'], headers=headers, cutout_datafolder=cutout_dir,
                               params={'limit': illustris_cfg['top_n'], 'primary_flag': 1,
                                       'order_by': '-mass_gas'})  # 'order_by':'-mass_dm'
            subs_sim_r += subs_sim_z_r['results']

        subs_r += subs_sim_r
        print(f"Number of subs loaded for {name}", len(subs_sim_r))

    print(f"Total subs loaded", len(subs_r))

    # Download all the cutouts
    cutouts_path = get_cutouts(subs=subs_r, headers=headers, cutout_datafolder=cutout_dir,
                               fail_on_error=env_cfg["fail_on_error"])
    with open(cutouts_path, "rb") as f:
        sc = pickle.load(f)

    for i, mode in enumerate(illustris_cfg["modes"].keys()):
        if mode not in ("proj", "slice"):
            raise ValueError(f"Expected either 'proj' or 'slice' at index {i} but got mode '{mode}'!")

    argument_list = []
    for res in sc:
        sub = res["sub"]
        cutout = res["cutout"]
        if isinstance(cutout, tuple):
            cutout = save_cutout(cutout, cutout_dir)
        # Generate the fits file
        # We can change the normal to rotate the image
        cutout_args = (cutout, dataset_dir, sub, illustris_cfg['modes'], cloudy_emissivity_root, illustris_cfg['emin'],
                       illustris_cfg['emax'], illustris_cfg['width'],
                       illustris_cfg['resolutions'],
                       illustris_cfg['redshift'], illustris_cfg['overwrite'], illustris_cfg['create_preview'])
        argument_list.append(cutout_args)

    if env_cfg["debug"]:
        for cutout_args in argument_list:
            cutout_to_xray_fits(*cutout_args)
    else:
        mp_run(cutout_to_xray_fits, argument_list, mp_conf=mp_cfg)
