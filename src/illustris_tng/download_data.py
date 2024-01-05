import re
from pathlib import Path
from typing import Iterator, Optional, Tuple, Union

import requests
from requests import Response
from loguru import logger

_cutout_request = {
    "gas": "Coordinates,Density,ElectronAbundance,GFM_Metallicity,InternalEnergy,Masses,Velocities"
}

# Matches urls like 'http://www.tng-project.org/api/TNG300-1/snapshots/99/subhalos/538308/cutout.hdf5'
_regex = re.compile(
    r"http://www.tng-project.org/api/.*/snapshots/\d*/\w*/\d*/cutout.hdf5"
)


def _handle_cutout_name(cutout_url: str) -> Tuple[str, str, str]:
    # Create a unique filename, use the url as base if possible
    if _regex.match(cutout_url):
        parts = cutout_url.split("/")[-6:-1]
        return parts[0], parts[2], f"{parts[3]}_{parts[4]}"
    else:
        raise ValueError


def get_available_simulations(api_key: str) -> Iterator[Tuple[str, str]]:
    content = get("https://www.tng-project.org/api/", headers={"api-key": api_key})
    simulations = content["simulations"]

    return (
        (simulation.pop("name"), simulation.pop("url")) for simulation in simulations
    )


def get_subhalos(
    api_key: str, simulation_url: str, snapshot_num: int, params: dict
) -> Iterator[str]:
    content = get(
        f"{simulation_url}/snapshots/{snapshot_num}/subhalos/",
        params=params,
        headers={"api-key": api_key},
    )
    subhalos = content["results"]

    return (subhalo.pop("url") for subhalo in subhalos)


def get(path: str, headers, params=None) -> Union[dict, Response]:
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers["content-type"] == "application/json":
        return r.json()  # parse json responses automatically

    return r


def get_cutouts(
    subhalo_url: str,
    api_key: str,
    cutout_datafolder: Path,
    fail_on_error: bool = False,
    cutouts_compressed: Optional[Path] = None,
) -> Optional[dict]:
    try:
        subhalo = get(subhalo_url, headers={"api-key": api_key})
        cutout_url = f"{subhalo_url}cutout.hdf5"
        filename = _handle_cutout_name(cutout_url=cutout_url)
        cutout_datafolder = cutout_datafolder / filename[0] / filename[1]
        cutout_datafolder.mkdir(parents=True, exist_ok=True)
        filename = filename[2]
        cutout_file = cutout_datafolder / f"{filename}.hdf5"

        if not cutout_file.exists():
            with requests.get(
                cutout_url,
                headers={"api-key": api_key},
                params=_cutout_request,
                stream=True,
            ) as r:
                r.raise_for_status()
                with open(cutout_file, "wb") as f:
                    for chunk in r.iter_content(chunk_size=int(1e6)):
                        f.write(chunk)

        cutout_dict = {
            "file": cutout_file.resolve(),
            "x": subhalo["pos_x"],
            "y": subhalo["pos_y"],
            "z": subhalo["pos_z"],
        }

        return cutout_dict

    except Exception as e:
        if fail_on_error:
            logger.exception(f"Failed to load sub {subhalo_url}!")
            raise
        else:
            logger.warning(f"Failed to load sub {subhalo_url} due to error {e}")
