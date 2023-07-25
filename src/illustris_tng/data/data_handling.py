import datetime
import json
import re
from pathlib import Path
from typing import Optional, Dict, Tuple

_reg_filename = "reg.json"
# Matches urls like 'http://www.tng-project.org/api/TNG300-1/snapshots/99/subhalos/538308/cutout.hdf5'
_regex = re.compile(r"http://www.tng-project.org/api/.*/snapshots/\d*/\w*/\d*/cutout.hdf5")


def write_dict_to_file(d: dict, path: Path) -> None:
    with open(path, "w") as file:
        json.dump(d, file)


def load_dict_from_file(path: Path) -> dict:
    with open(path, "r") as file:
        return json.load(file)


def update_reg_file(reg_path: Path, filename: str, cutout_url: str) -> None:
    filename_dict = load_dict_from_file(reg_path)
    filename_dict[cutout_url] = filename
    write_dict_to_file(filename_dict, reg_path)


# Since the cutout files are quite large (+- 1gb) we do not want to download them if we already have them
# in order to do this we check

def get_saved_file(cutout_url: str, datafolder: Path) -> Optional[str]:
    reg_path = datafolder / _reg_filename
    datafolder.mkdir(parents=True, exist_ok=True)
    filename_dict: Dict[str, str] = {}  # Start with an empty dict

    if reg_path.exists():
        # If the reg file exists, then load it
        filename_dict = load_dict_from_file(reg_path)

    return filename_dict.get(cutout_url, None)


def handle_cutout_name(cutout_url: str) -> str:
    # Create a unique filename, use the url as base if possible
    if _regex.match(cutout_url):
        parts = cutout_url.split('/')[-6:-1]
        return f"{parts[0]}_z_{parts[2]}_{parts[3]}_{parts[4]}.hdf5"  # Skips "snapshot"
    else:
        # The url is not the right format, use the datetime timestamp as name
        d = datetime.datetime.utcnow()
        epoch = datetime.datetime(1970, 1, 1)
        t = f"{(d - epoch).total_seconds()}"
        return f"{t.replace('.', '')}.hdf5"


# def save_cutout(cutout_url: str, content, datafolder: Path) -> str:
#     reg_path = datafolder / _reg_filename
#     # Create a unique filename, use the url as base if possible
#     if _regex.match(cutout_url):
#         parts = cutout_url.split('/')[-6:-1]
#         filename = f"{parts[0]}_z_{parts[2]}_{parts[3]}_{parts[4]}.hdf5"  # Skips "snapshot"
#     else:
#         # The url is not the right format, use the datetime timestamp as name
#         d = datetime.datetime.utcnow()
#         epoch = datetime.datetime(1970, 1, 1)
#         t = f"{(d - epoch).total_seconds()}"
#         filename = f"{t.replace('.', '')}.hdf5"
#
#     # Write the file
#     with open(datafolder / filename, "wb") as f:
#         f.write(content)
#
#     if reg_path.exists():
#         # Update the reg dict
#         update_reg_file(reg_path, filename, cutout_url)
#     else:
#         write_dict_to_file({cutout_url: filename}, reg_path)
#
#     return filename

def save_cutout(cutout: Tuple[Path, str, bytes], datafolder: Path) -> Path:
    filename, path, content = cutout[0], cutout[1], cutout[2]
    reg_path = datafolder / _reg_filename
    with open(datafolder / filename, "wb") as f:
        f.write(content)
    if reg_path.exists():
        # Update the reg dict
        update_reg_file(reg_path, str(filename.resolve()), path)
    else:
        write_dict_to_file({path: str(filename.resolve())}, reg_path)
    return datafolder / filename
