import tomllib
from argparse import ArgumentParser
from datetime import datetime, timedelta
from pathlib import Path

from loguru import logger

from src.config import EnergyCfg, EnvironmentCfg, MultiprocessingCfg, SimputCfg
from src.xmm_utils.run_utils import configure_logger, load_satellites
from src.xmm_utils.tools import generate_simput

logger.remove()

if __name__ == "__main__":
    parser = ArgumentParser(prog="", description="")
    parser.add_argument(
        "-a",
        "--agn_counts_file",
        default=Path(__file__).parent.resolve() / "res" / "agn_counts.cgi",
        type=Path,
        help="Path to agn_counts_cgi.",
    )
    parser.add_argument(
        "-p",
        "--config_path",
        type=Path,
        default=Path(__file__).parent.resolve() / "config.toml",
        help="Path to config file.",
    )

    args = parser.parse_args()

    with open(args.config_path, "rb") as file:
        cfg: dict[str, dict] = tomllib.load(file)
    env_cfg = EnvironmentCfg(**cfg.pop("environment"))

    simput_cfg = SimputCfg(
        **cfg.pop("simput"),
        simput_dir=env_cfg.working_dir / "simput",
        fits_dir=env_cfg.working_dir / "fits",
        fits_compressed=env_cfg.output_dir / "fits.tar.gz",
    )
    energies = EnergyCfg(**cfg.pop("energy"))
    mp_cfg = MultiprocessingCfg(**cfg.pop("multiprocessing"))

    satellites = load_satellites(cfg.pop("instruments"))

    del cfg

    configure_logger(
        log_dir=env_cfg.log_dir,
        log_name="02_generate_simput.log",
        enqueue=True,
        debug=env_cfg.debug,
        verbose=env_cfg.verbose,
        rotation=timedelta(hours=1),
        retention=2,
    )

    logger.info(f"Found satellites with instruments: {satellites}")

    starttime = datetime.now()
    generate_simput(
        simput_cfg=simput_cfg,
        energies=energies,
        env_cfg=env_cfg,
        mp_cfg=mp_cfg,
        satellites=satellites,
        agn_counts_file=args.agn_counts_file,
        delete_product=True,
    )
    endtime = datetime.now()
    logger.info(f"Duration: {endtime - starttime}")
