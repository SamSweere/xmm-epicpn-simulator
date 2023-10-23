from argparse import ArgumentParser
from pathlib import Path

from loguru import logger

from src.illustris_tng import run

logger.remove()

if __name__ == '__main__':
    parser = ArgumentParser(prog="", description="")
    parser.add_argument("-k", "--api_key", type=str, required=True,
                        help="IllustrisTNG API key. If you don't have one, create an account at "
                             "https://www.tng-project.org/data/")
    parser.add_argument("-p", "--config_path", type=Path, required=True, help="Path to config file.")
    parser.add_argument("-c", "--cloudy_emissivity", type=Path, required=True,
                        help="Path to root directory of cloudy_emissivity_v2.h5")

    args = parser.parse_args()
    run(path_to_cfg=args.config_path, api_key=args.api_key, cloudy_emissivity_root=args.cloudy_emissivity)
