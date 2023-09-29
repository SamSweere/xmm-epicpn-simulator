import json
from argparse import ArgumentParser
from pathlib import Path
from typing import Dict, Optional
from itertools import repeat
import numpy as np

from src.simput.gen import simput_generate
from src.xmm_utils.multiprocessing import mp_run


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

    working_directory = Path(env_cfg["working_directory"]).expanduser()

    tmp_dir = working_directory / "simput_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    debug = env_cfg["debug"]
    keep_files = env_cfg["keep_files"]
    verbose = env_cfg["verbose"]

    simput_dir = working_directory / "simput"
    simput_dir.mkdir(parents=True, exist_ok=True)
    sim_in_dataset_dir = working_directory / simput_cfg["simput_in_image_dataset"]

    # run_dir = working_directory / "tmp"
    # run_dir.mkdir(parents=True, exist_ok=True)

    sample_num = simput_cfg["num_img_sample"]
    zoom_range = simput_cfg["zoom_img_range"]
    sigma_b_range = simput_cfg['sigma_b_img_range']
    offset_std = simput_cfg['offset_std']

    argument_list = []
    mode_dict: Dict[str, int] = simput_cfg["mode"]
    for mode, num in mode_dict.items():
        if num < -1:
            raise ValueError("num has to be >= -1")
        if num == 0:
            continue

        mode_dir = simput_dir / mode
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
                    continue

                zoom = np.round(np.random.uniform(low=zoom_range[0], high=zoom_range[1], size=to_generate), 2)
                sigma_b = np.round(np.random.uniform(low=sigma_b_range[0], high=sigma_b_range[1], size=to_generate), 2)
                offset_x = np.round(rng.normal(loc=-offset_std, scale=offset_std, size=to_generate), 2)
                offset_y = np.round(rng.normal(loc=-offset_std, scale=offset_std, size=to_generate), 2)

                img_settings = {
                    "img_path": in_file.resolve(),
                    "zoom": zoom,
                    "sigma_b": sigma_b,
                    "offset_x": offset_x,
                    "offset_y": offset_y
                }

                argument_list.append((mode, img_settings, tmp_dir, mode_dir, keep_files, verbose))
        elif mode == "agn":
            if debug:
                img_settings = {
                    "num": num,
                    "agn_counts_file": agn_counts_file
                }
                argument_list.append((mode, img_settings, tmp_dir, mode_dir, keep_files, verbose))
            else:
                img_settings = {
                    "num": 1,
                    "agn_counts_file": agn_counts_file
                }
                argument_list.extend(repeat((mode, img_settings, tmp_dir, mode_dir, keep_files, verbose), num))
        elif mode == "background" or mode == "exposure_map":
            if debug:
                img_settings = {
                    "num": num,
                    "spectrum_file": spectrum_file
                }
                argument_list.append((mode, img_settings, tmp_dir, mode_dir, keep_files, verbose))
            else:
                img_settings = {
                    "num": 1,
                    "spectrum_file": spectrum_file
                }
                argument_list.extend(repeat((mode, img_settings, tmp_dir, mode_dir, keep_files, verbose), num))
        else:
            if debug:
                argument_list.append((mode, num, tmp_dir, mode_dir, keep_files, verbose))
            else:
                argument_list.extend(repeat((mode, 1, tmp_dir, mode_dir, keep_files, verbose), num))

    if debug:
        for args_tuple in argument_list:
            simput_generate(*args_tuple)
    else:
        mp_run(simput_generate, argument_list, mp_cfg)


if __name__ == '__main__':
    parser = ArgumentParser(prog="", description="")
    parser.add_argument("-a", "--agn_counts_file", type=Path, help="Path to agn_counts_cgi.")
    parser.add_argument("-p", "--config_path", type=Path, required=True, help="Path to config file.")
    parser.add_argument("-s", "--spectrum_file", type=Path, help="Path to spectrum file.")

    args = parser.parse_args()
    run(path_to_cfg=args.config_path, agn_counts_file=args.agn_counts_file, spectrum_file=args.spectrum_file)
