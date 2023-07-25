import json
from pathlib import Path
from typing import List, Dict, Union, Any

import numpy as np

from illustris_tng.data import get, get_cutouts
from illustris_tng.data.data_handling import save_cutout
from illustris_tng.gen.fits import cutout_to_xray_fits
from utils.multiprocessing import run_apply_async_multiprocessing


# Change the run directory to illustris_tng such that the illustris code finde the cloudy_emissivity_v2.h5 file
# root_dir = os.path.dirname(__file__)
# run_dir = tng_api_key_path = os.path.join(root_dir, 'illustris_tng')
# os.chdir(run_dir)
def run(path_to_cfg: Path):
    with open(path_to_cfg, "r") as file:
        cfg: Dict[str, dict] = json.load(file)
    env_cfg: Dict[str, Any] = cfg["environment"]
    mp_cfg: Dict[str, Any] = cfg["multiprocessing"]
    illustris_cfg: Dict[str, Any] = cfg["illustris"]

    api_key = illustris_cfg.get("api_key", None)

    if api_key is None:
        raise ValueError(f"API key is not set! If you do not have one, request it at "
                         f"https://www.tng-project.org/users/register/")

    working_dir = Path(env_cfg["working_directory"]).expanduser()
    cutout_dir = working_dir / "cutout_data"
    dataset_dir = working_dir / illustris_cfg["dataset_dir"]
    dataset_dir.mkdir(parents=True, exist_ok=True)

    baseUrl = "http://www.tng-project.org/api/"
    headers = {"api-key": api_key}

    print("Getting snaps...")
    simulations: List[Dict[str, Union[str, int]]] = get(baseUrl, headers=headers, cutout_datafolder=cutout_dir)[
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
    # TODO This depends on the amount of CPUs available. Maybe do this in run_apply_async_multiprocessing?
    chunks = np.array_split(subs_r, 12)
    argument_list = [(list(chunk), headers, cutout_dir) for chunk in chunks]
    sc = run_apply_async_multiprocessing(get_cutouts, argument_list, 12)
    # sc = get_cutouts(subs=subs_r, headers=headers, cutout_datafolder=cutout_dir)
    # sc = itertools.chain.from_iterable(sc)
    argument_list = []

    for res in sc:
        sub = res["sub"]
        cutout = res["cutout"]
        if isinstance(cutout, tuple):
            cutout = save_cutout(cutout, cutout_dir)
        # Generate the fits file
        # We can change the normal to rotate the image
        for mode in illustris_cfg['modes']:
            if mode == "proj":
                normals = illustris_cfg['proj_normals']
            elif mode == "slice":
                normals = illustris_cfg['slice_axes']
            else:
                raise ValueError(f"mode {mode} not suported")

            for width in illustris_cfg['width']:
                for normal in normals:
                    for resolution in illustris_cfg['resolutions']:
                        args = (
                            cutout, dataset_dir, sub, mode, illustris_cfg['emin'], illustris_cfg['emax'],
                            normal, width, resolution, illustris_cfg['redshift'], illustris_cfg['overwrite'],
                            illustris_cfg['create_preview'])
                        argument_list.append(args)

            # For debugging, this does not use multiprocessing
            # cutout_to_xray_fits(cutout, cutout_datafolder, dataset_dir, sub, mode, emin, emax,
            #                                             normal, width, resolution, redshift, overwrite,
            #                                             create_preview)

    run_apply_async_multiprocessing(cutout_to_xray_fits, argument_list, num_processes=12,  # TODO
                                    gb_per_process=mp_cfg['gb_per_process'])


if __name__ == '__main__':
    run(Path("/home/bojantodorkov/Projects/xmm-epicpn-simulator/cfg/illustris.json"))
