import pathlib
import tomllib
from argparse import ArgumentParser
from pathlib import Path

from loguru import logger

from src.config import EnergyCfg, EnvironmentCfg, MultiprocessingCfg, SimulationCfg
from src.tools.run_utils import configure_logger, load_satellites
from src.tools.tools import run_simulations

logger.remove()


if __name__ == "__main__":
    parser = ArgumentParser(prog="", description="")
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
    sim_cfg = SimulationCfg(
        **cfg.pop("simulation"),
        simput_dir=env_cfg.working_dir / "simput",
        out_dir=env_cfg.working_dir / "xmm_sim_dataset",
    )
    energies = EnergyCfg(**cfg.pop("energy"))
    mp_cfg = MultiprocessingCfg(**cfg.pop("multiprocessing"))

    satellites = load_satellites(cfg.pop("instruments"))

    del cfg

    configure_logger(
        log_dir=env_cfg.log_dir,
        log_name="03_xmm_simulation.log",
        enqueue=True,
        debug=env_cfg.debug,
        verbose=env_cfg.verbose,
    )

    run_simulations(
        sim_cfg=sim_cfg,
        energies=energies,
        env_cfg=env_cfg,
        mp_cfg=mp_cfg,
        satellites=satellites,
    )
