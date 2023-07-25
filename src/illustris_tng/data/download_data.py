from pathlib import Path
from typing import Union, List, Tuple
from warnings import warn

import requests
from tqdm import tqdm

from illustris_tng.data.data_handling import get_saved_file, save_cutout, handle_cutout_name

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

    if 'cutout' in path and 'content-disposition' in r.headers:
        filename = save_cutout(cutout_url=path, content=r.content, datafolder=cutout_datafolder)
        # return filename  # return the filename string
        return path, r.content, cutout_datafolder
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


def get_cutouts(subs, headers, cutout_datafolder) -> List[dict]:
    sc = []
    for sub in tqdm(subs):
        try:
            sub = _get_sub(sub["url"], headers=headers)
            cutout = _get_cutout(sub["cutouts"]["subhalo"], headers=headers, cutout_datafolder=cutout_datafolder)
            sc.append({
                "sub": sub,
                "cutout": cutout
            })
        except Exception as e:
            warn(f"Failed to load sub {sub['url']} due to error {e}")
    return sc
