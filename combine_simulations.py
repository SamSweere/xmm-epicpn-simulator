# This file should not depend on other in project dependencies
import json
import os
import random
from datetime import datetime
from itertools import repeat
from pathlib import Path
from typing import List, Union, Dict, Optional, Tuple

import numpy as np
from astropy.io import fits

config = {
    'run_dir': None,
    'data_dir': None,
    'dataset_dir': None,
    'dataset_name': 'xmm_demo_dataset',
}


def process_one_simulation(file_path: Path, general_header_dict: Optional[dict]) -> Tuple[np.ndarray, dict, dict]:
    with fits.open(file_path) as hdu:
        hdu_data = hdu['PRIMARY'].data
        header = hdu['PRIMARY'].header

        hdu_general_header_dict = {
            'TELESCOP': (header['TELESCOP'], header.comments['TELESCOP']),
            'WCSAXES': (header['WCSAXES'], header.comments['WCSAXES']),
            'CRPIX1': (header['CRPIX1'], header.comments['CRPIX1']),
            'CRPIX2': (header['CRPIX2'], header.comments['CRPIX2']),
            'CDELT1': (header['CDELT1'], header.comments['CDELT1']),
            'CDELT2': (header['CDELT2'], header.comments['CDELT2']),
            'CUNIT1': (header['CUNIT1'], header.comments['CUNIT1']),
            'CUNIT2': (header['CUNIT2'], header.comments['CUNIT2']),
            'CTYPE1': (header['CTYPE1'], header.comments['CTYPE1']),
            'CTYPE2': (header['CTYPE2'], header.comments['CTYPE2']),
            'CRVAL1': (header['CRVAL1'], header.comments['CRVAL1']),
            'CRVAL2': (header['CRVAL2'], header.comments['CRVAL2']),
            'LONPOLE': (header['LONPOLE'], header.comments['LONPOLE']),
            'RADESYS': (header['RADESYS'], header.comments['RADESYS']),
            'EXPOSURE': (header['EXPOSURE'], header.comments['EXPOSURE']),
            'RESMULT': (header['RESMULT'], header.comments['RESMULT'])
        }

    if general_header_dict is None:
        general_header_dict = hdu_general_header_dict
    else:
        if general_header_dict != hdu_general_header_dict:
            raise ValueError(f"General header values are not the same, something is wrong.\n"
                             f"general_header_dict = {general_header_dict}\n"
                             f"hdu_general_header_dict = {hdu_general_header_dict}")

    return hdu_data, header, general_header_dict


def apply_det_mask(data, det_mask_path) -> np.ndarray:
    with fits.open(det_mask_path) as det_mask:
        data = data * det_mask[0].data

    return data


def combine_simulation(
        img_path: Path,
        agn_path: Path = None,
        background_path: Path = None,
        det_mask_path: Path = None
) -> Tuple[np.ndarray, fits.Header]:
    if not img_path and not agn_path and not background_path:
        raise AssertionError("Not all inputfiles can be None")

    header = fits.Header()

    data, img_header, general_header_dict = process_one_simulation(img_path, None)

    header["IMG_FILE"] = (img_path.name, "Input simulation file used as img")
    header["IMG_SIMP"] = (img_header['SIMPUT'], "Input img simput used for simulation")

    if agn_path is not None:
        agn_data, agn_header, general_header_dict = process_one_simulation(agn_path, general_header_dict)

        header["AGN_FILE"] = (agn_path.name, "Input simulation file used as agn")
        header["AGN_SIMP"] = (agn_header['SIMPUT'], "Input agn simput used for simulation")

        data = data + agn_data

    if background_path is not None:
        background_data, background_header, general_header_dict = process_one_simulation(background_path,
                                                                                         general_header_dict)

        header["BCK_FILE"] = (background_path.name, "Input simulation file used as background")
        header["BCK_SIMP"] = (background_header['SIMPUT'], "Input background simput used for simulation")

        data = data + background_data

    if det_mask_path is not None:
        # Apply the detector mask
        data = apply_det_mask(data, det_mask_path)
        header["DETMASK"] = (det_mask_path.name, "Detector mask used on this image")

    # Add the general background to the header
    header.update(general_header_dict)

    header['COMMENT'] = f"Created by Sam Sweere (samsweere@gmail.com) for ESAC at " \
                        f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}"

    return data, header


def random_sample_sim_files(sample_dir: Path, num: int) -> List[Path]:
    files = list(sample_dir.glob("*.fits.gz"))
    files.extend(sample_dir.glob("*.fits"))

    if len(files) < num:
        print(f"WARNING: sampling {num} files from a source of {len(files)} files. This will cause (many) duplicates")

    sampled_files = random.sample(files, num)

    return sampled_files


def combine_and_save_sim(
        output_dir: Path,
        imgs: List[Path],
        agns: List[Path] = None,
        backgrounds: List[Path] = None,
        det_mask: Path = None
):
    if not os.path.exists(output_dir):
        raise NotADirectoryError(f"Output directory {output_dir} does not exist")

    if agns is None:
        agns = repeat(None, len(imgs))

    if backgrounds is None:
        backgrounds = repeat(None, len(imgs))

    for img, agn, background in zip(imgs, agns, backgrounds):
        print(f"Combining: img={img}, agn={agn} and background={background}")
        data, header = combine_simulation(img_path=img, agn_path=agn, background_path=background,
                                          det_mask_path=det_mask)

        # Determine the output name
        out_name = img.name.replace(".fits", "").replace(".gz", "")
        if agn is not None:
            agn_name = agn.name.replace(".fits", "").replace(".gz", "")
            agn_name = agn_name.split("_")[1]
            out_name = f"{out_name}_agn_{agn_name}"

        if background is not None:
            background_name = background.name.replace(".fits", "").replace(".gz", "")
            background_name = background_name.split("_")[1]
            out_name = f"{out_name}_back_{background_name}"

        compressed_out_path = output_dir / f"{out_name}.fits.gz"

        hdu = fits.PrimaryHDU(data=data, header=header)
        hdu.writeto(compressed_out_path, overwrite=True)


def run(cfg: Union[Path, Dict[str, dict]]) -> None:
    if isinstance(cfg, Path):
        with open(cfg, "r") as f:
            cfg: Dict[str, dict] = json.load(f)

    env_cfg = cfg["environment"]

    exposures = cfg["exposure"]
    mode_dict = cfg["mode"]
    res_mults = cfg['res_mult']

    lowest_res_mult = f"{res_mults[0]}x"
    working_directory = Path(env_cfg["working_directory"]).expanduser()
    dataset_dir = working_directory / env_cfg["dataset_dir"]

    for exposure in exposures:
        exposure = f"{exposure}ks"
        print(f"Processing exposure: {exposure}")
        exposure_dir = dataset_dir / exposure
        for mode, mode_params in mode_dict.items():
            amount = mode_params["amount"]
            if mode == "background":
                amount_agn = 0
                amount_bg = 0
            else:
                amount_agn = mode_params["amount_agn"]
                amount_bg = mode_params["amount_bg"]

            if amount < -1:
                raise ValueError(f"amount of {mode} has to be an int >= -1, but got {amount}!")
            if amount_agn < 0:
                raise ValueError(f"amount_agn of {mode} has to be an int >= 0, but got {amount_agn}")
            if amount_bg < 0:
                raise ValueError(f"amount_bg of {mode} has to be an int >= 0, but got {amount_bg}!")

            if amount == 0:
                continue

            print(f"Processing mode: {mode}")

            for res_mult in res_mults:
                res_mult = f"{res_mult}x"
                print(f"Processing res mult: {res_mult}")

                in_dir = exposure_dir / mode / res_mult
                output_dir = dataset_dir / "combined" / exposure / mode / res_mult
                output_dir.mkdir(exist_ok=True, parents=True)

                imgs = list(in_dir.glob("*.fits.gz"))
                imgs.extend(in_dir.glob("*.fits"))
                if amount != -1:
                    imgs = imgs[:amount]

                # Get the detector mask path, these should be present in the simulated dataset
                det_mask = None
                if mode_params["add_detmask"]:
                    # TODO Fix path for different kinds of detectors
                    det_mask = Path.cwd() / "res" / "detector_masks" / f"pn_mask_500_2000_detxy_{res_mult}.ds"

                agn_files = None
                if amount_agn:
                    agn_files = random_sample_sim_files(sample_dir=exposure_dir / "agn" / res_mult, num=amount_agn)

                background_files = None
                if amount_bg:
                    background_files = random_sample_sim_files(sample_dir=exposure_dir / "background" / res_mult,
                                                               num=amount_bg)

                combine_and_save_sim(output_dir=output_dir, imgs=imgs, agns=agn_files, backgrounds=background_files,
                                     det_mask=det_mask)


if __name__ == '__main__':
    run(Path("/home/bojantodorkov/Projects/xmm-epicpn-simulator/cfg/combine.json"))
