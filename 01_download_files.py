import os
import pathlib
import tomllib
from argparse import ArgumentParser
from datetime import datetime, timedelta
from pathlib import Path

from dotenv import load_dotenv
from loguru import logger

from src.config import DownloadCfg, EnergyCfg, EnvironmentCfg, MultiprocessingCfg
from src.xmm_utils.run_utils import configure_logger
from src.xmm_utils.tools import download_data

logger.remove()

if __name__ == "__main__":
    # Load .env file
    load_dotenv()

    # Get TNG_API_KEY from environment variables, or None if it's not present
    default_tng_api_key = os.getenv("TNG_API_KEY")

    parser = ArgumentParser(prog="", description="")
    # Add TNG_API_KEY argument
    parser.add_argument(
        "-k",
        "--api_key",
        type=str,
        default=default_tng_api_key,
        required=default_tng_api_key is None,
        help="TNG API key. If you don't have one, create an account at " "https://www.tng-project.org/data/",
    )

    parser.add_argument(
        "-p",
        "--config_path",
        type=Path,
        default=pathlib.Path(__file__).parent.resolve() / "config.toml",
        help="Path to config file.",
    )

    args = parser.parse_args()

    with open(args.config_path, "rb") as file:
        cfg: dict[str, dict] = tomllib.load(file)
    env_cfg = EnvironmentCfg(**cfg.pop("environment"))

    download_cfg = DownloadCfg(
        **cfg.pop("download"),
        cutouts_path=env_cfg.working_dir / "cutouts",
        cutouts_compressed=env_cfg.output_dir / "cutouts.tar.gz",
        fits_path=env_cfg.working_dir / "fits",
        fits_compressed=env_cfg.output_dir / "fits.tar.gz",
    )
    energies = EnergyCfg(**cfg.pop("energy"))
    mp_cfg = MultiprocessingCfg(**cfg.pop("multiprocessing"))

    del cfg

    configure_logger(
        log_dir=env_cfg.log_dir,
        log_name="01_download_files.log",
        enqueue=True,
        debug=env_cfg.debug,
        verbose=env_cfg.verbose,
        rotation=timedelta(hours=1),
        retention=2,
    )

    starttime = datetime.now()
    download_data(
        download_cfg=download_cfg,
        energies=energies,
        env_cfg=env_cfg,
        mp_cfg=mp_cfg,
        api_key=args.api_key,
        delete_product=True,
    )
    endtime = datetime.now()
    logger.info(f"Duration: {endtime - starttime}")
