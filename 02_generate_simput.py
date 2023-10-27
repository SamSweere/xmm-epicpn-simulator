import json
from argparse import ArgumentParser
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Dict, Optional

import numpy as np
from loguru import logger

from src.simput import constants
from src.simput.gen import simput_generate
from src.xmm_utils.multiprocessing import get_num_processes

logger.remove()


def handle_error(error):
    logger.exception(error)


def run(
        path_to_cfg: Path,
        agn_counts_file: Optional[Path],
        spectrum_file: Optional[Path]
) -> None:
    with open(path_to_cfg, "r") as f:
        cfg: Dict[str, dict] = json.load(f)
    env_cfg = cfg["environment"]
    mp_cfg = cfg["multiprocessing"]
    simput_cfg = cfg["simput"]

    log_dir = Path(env_cfg["log_directory"]).expanduser()
    log_dir.mkdir(parents=True, exist_ok=True)
    log_level = "DEBUG" if env_cfg["debug"] else "INFO"
    log_file = log_dir / "02_generate_simput.log"
    logger.add(f"{log_file.resolve()}", enqueue=True,
               format="{time:YYYY-MM-DD at HH:mm:ss} | {level} | {message}", level=log_level)
    log_file.chmod(0o777)

    working_directory = Path(env_cfg["working_directory"]).expanduser()

    tmp_dir = working_directory / "simput_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    debug = env_cfg["debug"]
    verbose = env_cfg["verbose"]

    if not debug:
        logger.info(f"Since 'debug' is set to 'false' the generation will be run asynchronously.")

    simput_dir = working_directory / "simput"
    simput_dir.mkdir(parents=True, exist_ok=True)
    sim_in_dataset_dir = working_directory / simput_cfg["simput_in_image_dataset"]

    # run_dir = working_directory / "tmp"
    # run_dir.mkdir(parents=True, exist_ok=True)

    emin = simput_cfg["emin"]
    emax = simput_cfg["emax"]

    sample_num = simput_cfg["num_img_sample"]
    zoom_range = simput_cfg["zoom_img_range"]
    sigma_b_range = simput_cfg['sigma_b_img_range']
    offset_std = simput_cfg['offset_std']
    instrument_name = simput_cfg["instrument_name"]

    mode_dict: Dict[str, int] = simput_cfg["mode"]
    with Pool(get_num_processes(mp_conf=mp_cfg)) as pool:
        for mode, num in mode_dict.items():
            if mode not in constants.available_modes:
                raise ValueError(f"Unkown mode '{mode}'! Available modes: {constants.available_modes}.")

            if num < -1:
                raise ValueError("num has to be >= -1")

            if num == 0:
                continue

            mode_dir = simput_dir / instrument_name / mode
            mode_dir.mkdir(parents=True, exist_ok=True)

            if mode == "img":
                in_files = list(sim_in_dataset_dir.glob("*.fits"))
                if not num == -1:
                    in_files = in_files[:num]

                rng = np.random.default_rng()

                for in_file in in_files:
                    generated_files = len(list(mode_dir.glob(f"{in_file.stem}*")))
                    to_generate = sample_num - generated_files

                    if to_generate <= 0:
                        logger.info(f"Won't generate any images for {in_file.name}.")
                        continue

                    logger.info(f"Will generate {to_generate} images for {in_file.name}.")

                    zoom = np.round(rng.uniform(low=zoom_range[0], high=zoom_range[1], size=to_generate), 2)
                    sigma_b = np.round(rng.uniform(low=sigma_b_range[0], high=sigma_b_range[1], size=to_generate), 2)
                    offset_x = np.round(rng.normal(loc=-offset_std, scale=offset_std, size=to_generate), 2)
                    offset_y = np.round(rng.normal(loc=-offset_std, scale=offset_std, size=to_generate), 2)

                    img_settings = {
                        "img_path": in_file.resolve(),
                        "zoom": zoom,
                        "sigma_b": sigma_b,
                        "offset_x": offset_x,
                        "offset_y": offset_y
                    }

                    arguments = (instrument_name, emin, emax, mode, img_settings, tmp_dir, mode_dir, verbose)
                    if debug:
                        pool.apply(simput_generate, arguments)
                    else:
                        pool.apply_async(simput_generate, arguments, error_callback=handle_error)

            else:
                if mode == "agn":
                    logger.info(f"Will generate {num} AGNs.")
                    img_settings = {"agn_counts_file": agn_counts_file}
                elif mode == "background" or mode == "exposure_map":
                    logger.info(f"Will generate {num} {'backgrounds' if mode == 'background' else 'exposure maps'}.")
                    img_settings = {"spectrum_file": spectrum_file}
                else:
                    logger.info(f"Will generate {num} {'test grids' if mode == 'test_grid' else 'random sources'}.")
                    img_settings = {}

                if debug:
                    img_settings["num"] = num
                    arguments = (instrument_name, emin, emax, mode, img_settings, tmp_dir, mode_dir, verbose)
                    pool.apply(simput_generate, arguments)
                else:
                    img_settings["num"] = 1
                    for _ in range(num):
                        arguments = (instrument_name, emin, emax, mode, img_settings, tmp_dir, mode_dir, verbose)
                        pool.apply_async(simput_generate, arguments, error_callback=handle_error)
        pool.close()
        pool.join()
    logger.info("Done!")


if __name__ == '__main__':
    parser = ArgumentParser(prog="", description="")
    parser.add_argument("-a", "--agn_counts_file", type=Path, help="Path to agn_counts_cgi.")
    parser.add_argument("-p", "--config_path", type=Path, required=True, help="Path to config file.")
    parser.add_argument("-s", "--spectrum_file", type=Path, help="Path to spectrum file.")

    args = parser.parse_args()
    run(path_to_cfg=args.config_path, agn_counts_file=args.agn_counts_file, spectrum_file=args.spectrum_file)
