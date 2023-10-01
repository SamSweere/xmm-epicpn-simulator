import json
from argparse import ArgumentParser
from pathlib import Path
from typing import Dict, Union

from src.sixte.simulator import run_simulation_modes


def run(path_to_cfg: Union[Path, Dict[str, dict]]) -> None:
    if isinstance(path_to_cfg, Path):
        with open(path_to_cfg, "r") as f:
            cfg: Dict[str, dict] = json.load(f)

    env_cfg = cfg["environment"]
    mp_cfg = cfg["multiprocessing"]
    instrument_cfg = cfg["instrument"]

    working_directory = Path(env_cfg["working_directory"]).expanduser()

    dataset_dir = working_directory / "xmm_sim_dataset"
    dataset_dir.mkdir(parents=True, exist_ok=True)

    simput_path = working_directory / "simput"
    if not simput_path.exists():
        raise NotADirectoryError(f"Simput directory '{simput_path.resolve()}' not found! "
                                 f"Please make sure that the simput files have been generated and that "
                                 f"the working directory is set correctly.")

    mode_amount_dict: Dict[str, int] = instrument_cfg["mode"]
    run_simulation_modes(mp_conf=mp_cfg,
                         mode_amount_dict=mode_amount_dict,
                         working_directory=working_directory,
                         simput_path=simput_path,
                         instrument_conf=instrument_cfg,
                         debug=env_cfg["debug"],
                         verbose=env_cfg["verbose"])


if __name__ == '__main__':
    parser = ArgumentParser(prog="", description="")
    parser.add_argument("-p", "--config_path", type=Path, required=True, help="Path to config file.")

    args = parser.parse_args()
    run(path_to_cfg=args.config_path)
