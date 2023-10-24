from pathlib import Path
from typing import Union, List, Tuple

import requests
from loguru import logger

from src.illustris_tng.data.data_handling import get_saved_file, handle_cutout_name
from src.illustris_tng.data.data_handling import save_cutout

_cutout_request = {'gas': 'Coordinates,Density,ElectronAbundance,GFM_Metallicity,'
                          'InternalEnergy,Masses,NeutralHydrogenAbundance,Velocities'}


def get(
        path: str,
        headers,
        cutout_datafolder: Path,
        params=None
) -> Union[str, requests.Response, dict, List[dict]]:
    if 'cutout' in path:
        # If the path is a cutout getter, check if we already have it
        filename = get_saved_file(path, cutout_datafolder)
        if filename:
            # We have the file return the name
            return filename

    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json()  # parse json responses automatically

    return r


def _get_sub(path, headers):
    response = requests.get(path, headers=headers)
    response.raise_for_status()
    return response.json()


def _get_cutout(path, headers, cutout_datafolder: Path) -> Union[Path, Tuple[Path, str, bytes]]:
    filename = get_saved_file(path, cutout_datafolder)
    if filename is not None:
        return cutout_datafolder / filename
    else:
        r = requests.get(path, headers=headers, params=_cutout_request)
        r.raise_for_status()
        filename = handle_cutout_name(path)
        return cutout_datafolder / filename, path, r.content


def get_cutouts(
        subs,
        headers,
        cutout_datafolder: Path,
        fail_on_error: bool = False
) -> List[dict]:
    sc = []
    for sub in subs:
        url = sub["url"]
        try:
            sub = _get_sub(url, headers=headers)
            cutout = _get_cutout(sub["cutouts"]["subhalo"], headers=headers, cutout_datafolder=cutout_datafolder)
            if isinstance(cutout, tuple):
                cutout = save_cutout(cutout=cutout, datafolder=cutout_datafolder)
            sub = {
                "pos_x": sub["pos_x"],
                "pos_y": sub["pos_y"],
                "pos_z": sub["pos_z"]
            }
            sc.append({
                "sub": sub,
                "cutout": cutout
            })
        except Exception as e:
            if fail_on_error:
                logger.exception(f"Failed to load sub {url}!")
            else:
                logger.warning(f"Failed to load sub {url} due to error {e}")
    return sc
