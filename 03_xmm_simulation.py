import os
import pathlib
import shutil
import tomllib
from argparse import ArgumentParser
from datetime import datetime, timedelta
from functools import partial
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal

from loguru import logger

from src.config import EnergySettings, EnvironmentCfg, SimulationCfg
from src.sixte.simulator import run_xmm_simulation
from src.xmm.utils import create_mask, create_psf_file, create_vinget_file, create_xml_files
from src.xmm_utils.file_utils import compress_targz, decompress_targz
from src.xmm_utils.multiprocessing import mp_run
from src.xmm_utils.run_utils import configure_logger, load_satellites

logger.remove()


def _simulate_mode(
    instrument_name: Literal["epn", "emos1", "emos2"],
    max_event_pattern: int,
    mode: str,
    amount: int,
    sim_dir: Path,
    xmm_filter: str,
    sim_separate_ccds: bool,
    xmm_filter_dir: Path,
    xml_dir: Path,
) -> None:
    logger.info(f"START\tSimulating {instrument_name} for {mode.upper()}.")

    # Find the simput files
    mode_dir = sim_cfg.simput_dir / mode
    if mode == "agn":
        mode_dir = mode_dir / instrument_name

    if mode != "bkg":
        mode_glob = mode_dir.rglob("*.simput.gz")
        simputs = []
        for i, simput in enumerate(mode_glob, start=1):
            simputs.append(simput)
            if amount != -1 and i >= amount:
                break
    else:
        simputs = [next(mode_dir.rglob(f"*{instrument_name}.simput.gz"))] * amount

    to_run = partial(
        run_xmm_simulation,
        instrument_name=instrument_name,
        xml_dir=xml_dir.resolve(),
        mode=f"{mode}",
        tmp_dir=sim_dir.resolve(),
        out_dir=xmm_filter_dir.resolve(),
        max_event_pattern=max_event_pattern,
        exposure=sim_cfg.max_exposure,
        xmm_filter=xmm_filter,
        sim_separate_ccds=sim_separate_ccds,
        consume_data=env_cfg.consume_data,
    )
    kwds = (
        {"simput_file": simput.resolve(), "res_mult": res_mult, "emask": emasks[instrument_name][res_mult]}
        for res_mult in sim_cfg.res_mults
        for simput in simputs
    )
    _, duration = mp_run(to_run, kwds, sim_cfg.num_processes, env_cfg.debug)
    logger.success(f"DONE\tSimulating {instrument_name} for {mode.upper()}. Duration: {duration}")

    if env_cfg.working_dir != env_cfg.output_dir:
        mode_compressed = env_cfg.output_dir / "xmm_sim_dataset" / instrument_name / xmm_filter / f"{mode}.tar.gz"
        mode_compressed.parent.mkdir(parents=True, exist_ok=True)

        logger.info(
            f"Simulated {mode.upper()} will be compressed and moved to {mode_compressed.resolve()}."
            "Existing file will be overwritten."
        )

        compress_targz(
            in_path=xmm_filter_dir / f"{mode}",
            out_file_path=sim_dir / f"{mode}.tar.gz",
            remove_files=True,
        )
        shutil.move(src=sim_dir / f"{mode}.tar.gz", dst=mode_compressed)

        for xmm_mode_dir in xmm_filter_dir.rglob(f"{mode}{os.sep}"):
            shutil.rmtree(xmm_mode_dir)


def run(path_to_cfg: Path) -> None:
    starttime = datetime.now()
    with open(path_to_cfg, "rb") as file:
        cfg: dict[str, dict] = tomllib.load(file)

    global env_cfg, sim_cfg

    env_cfg = EnvironmentCfg(**cfg.pop("environment"))
    sim_cfg = SimulationCfg(
        **cfg.pop("simulation"),
        simput_dir=env_cfg.working_dir / "simput",
        out_dir=env_cfg.working_dir / "xmm_sim_dataset",
    )
    energies = EnergySettings(**cfg.pop("energy"))

    satellites = load_satellites(cfg.pop("instruments"))

    del cfg

    configure_logger(
        log_dir=env_cfg.log_dir,
        log_name="03_xmm_simulation.log",
        enqueue=True,
        debug=env_cfg.debug,
        verbose=env_cfg.verbose,
        rotation=timedelta(hours=1),
        retention=2,
    )

    with TemporaryDirectory(prefix="xml_") as xml_dir, TemporaryDirectory(prefix="sim_") as sim_dir:
        xml_dir = Path(xml_dir)
        sim_dir = Path(sim_dir)

        # Create all needed directories
        for sat in satellites:
            for name, instrument in sat:
                if not instrument.use:
                    continue
                filter_dir: Path = xml_dir / name / instrument.filter
                for res_mult in sim_cfg.res_mults:
                    (filter_dir / f"{res_mult}x").mkdir(exist_ok=True, parents=True)

        # Decrompress SIMPUT files if needed
        if env_cfg.working_dir != env_cfg.output_dir:
            for mode, amount in sim_cfg.modes:
                if amount == 0:
                    logger.debug(f"Skipping {mode} since simulation ammount is set to 0.")
                    continue

                simput_dir = env_cfg.output_dir / "simput" / mode

                simput_compressed_files = [next(simput_dir.rglob("*.tar.gz"))]

                for simput_compressed in simput_compressed_files:
                    if simput_compressed.exists():
                        logger.info(f"Found compressed SIMPUT files in {simput_compressed.resolve()}. Decompressing...")
                        decompress_targz(
                            in_file_path=simput_compressed,
                            out_file_dir=sim_cfg.simput_dir / mode,
                            tar_options="--strip-components=1",
                        )
                        logger.success("DONE\tDecompressing SIMPUT files.")

        logger.info("START\tCreating all PSF files.")
        to_run = partial(create_psf_file, xml_dir=xml_dir)

        kwds = (
            {"instrument_name": name, "res_mult": res_mult}
            for sat in satellites
            for name, instrument in sat
            if instrument.use
            for res_mult in sim_cfg.res_mults
        )
        _, duration = mp_run(to_run, kwds, sim_cfg.num_processes, env_cfg.debug)
        logger.success(f"DONE\tPSF files have been created. Duration: {duration}")

        logger.info("START\tCreating all vignetting files.")
        to_run = partial(create_vinget_file, xml_dir=xml_dir)

        kwds = ({"instrument_name": name} for name, instrument in sat if instrument.use)
        _, duration = mp_run(to_run, kwds, sim_cfg.num_processes, env_cfg.debug)
        logger.success(f"DONE\tVignetting files have been created. Duration: {duration}")

        logger.info("START\tCreating all XML files.")
        to_run = partial(
            create_xml_files,
            xml_dir=xml_dir,
            wait_time=sim_cfg.wait_time,
        )

        kwds = (
            {
                "instrument_name": name,
                "xmm_filter": instrument.filter,
                "sim_separate_ccds": instrument.sim_separate_ccds,
                "res_mult": res_mult,
            }
            for sat in satellites
            for name, instrument in sat
            if instrument.use
            for res_mult in sim_cfg.res_mults
        )
        _, duration = mp_run(to_run, kwds, sim_cfg.num_processes, env_cfg.debug)
        logger.success(f"DONE\tXML files have been created. Duration: {duration}")

        logger.info("START\tCreating all EMASKs.")
        to_run = partial(
            create_mask,
            observation_id="0935190401",
            res_mults=sim_cfg.res_mults,
            energies=energies,
        )

        kwds = (
            {
                "instrument_name": name,
                "mask_level": instrument.mask_level,
            }
            for sat in satellites
            for name, instrument in sat
            if instrument.use
        )
        global emasks
        emasks, duration = mp_run(to_run, kwds, sim_cfg.num_processes, env_cfg.debug)
        logger.success(f"DONE\tEMASKs have been created. Duration: {duration}")

        emasks = {key: value for d in emasks for key, value in d.items()}

        for mode, amount in sim_cfg.modes:
            if amount == 0:
                logger.info(f"Skipping {mode.upper()} simulation since amount is 0.")
                continue

            for sat in satellites:
                for name, instrument in sat:
                    if not instrument.use:
                        continue
                    xmm_filter_dir = sim_cfg.out_dir / name / instrument.filter
                    xmm_filter_dir.mkdir(exist_ok=True, parents=True)

                    _simulate_mode(
                        instrument_name=name,
                        max_event_pattern=instrument.max_event_pattern,
                        mode=mode,
                        amount=amount,
                        sim_dir=sim_dir,
                        xmm_filter=instrument.filter,
                        sim_separate_ccds=instrument.sim_separate_ccds,
                        xmm_filter_dir=xmm_filter_dir,
                        xml_dir=xml_dir,
                    )

            if env_cfg.working_dir != env_cfg.output_dir:
                shutil.rmtree(sim_cfg.simput_dir / f"{mode}")

    endtime = datetime.now()
    logger.info(f"Duration: {endtime - starttime}")


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
    run(path_to_cfg=args.config_path)
