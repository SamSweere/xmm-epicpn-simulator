import json
import os
import shutil
from argparse import ArgumentParser
from datetime import datetime, timedelta
from functools import partial
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal

from loguru import logger

from src.config import EnvironmentCfg, SimulationCfg
from src.sixte.simulator import run_xmm_simulation
from src.xmm.utils import create_psf_file, create_vinget_file, create_xml_files
from src.xmm_utils.file_utils import compress_targz, decompress_targz
from src.xmm_utils.multiprocessing import mp_run
from src.xmm_utils.run_utils import configure_logger

logger.remove()


def _simulate_mode(
    instrument_name: Literal["epn", "emos1", "emos2"],
    mode: str,
    amount: int,
    sim_dir: Path,
    xmm_filter_dir: Path,
    xml_dir: Path,
):
    logger.info(f"START\tSimulating {instrument_name} for {mode.upper()}.")
    mode_dir = sim_cfg.simput_dir / mode
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
        exposure=sim_cfg.max_exposure,
        xmm_filter=sim_cfg.filter,
        sim_separate_ccds=sim_cfg.sim_separate_ccds,
        consume_data=env_cfg.consume_data,
    )
    kwds = (
        {"simput_file": simput.resolve(), "res_mult": res_mult} for res_mult in sim_cfg.res_mults for simput in simputs
    )
    _, duration = mp_run(to_run, kwds, sim_cfg.num_processes, env_cfg.debug)
    logger.success(f"DONE\tSimulating {instrument_name} for {mode.upper()}. Duration: {duration}")

    if env_cfg.working_dir != env_cfg.output_dir:
        mode_compressed = env_cfg.output_dir / "xmm_sim_dataset" / instrument_name / sim_cfg.filter / f"{mode}.tar.gz"
        mode_compressed.parent.mkdir(parents=True, exist_ok=True)

        logger.info(
            f"Simulated {mode.upper()} will be compressed and moved to {mode_compressed.resolve()}."
            "Existing file will be overwritten."
        )

        compress_targz(
            in_path=xmm_filter_dir,
            out_file_path=sim_dir / f"{mode}.tar.gz",
            remove_files=True,
        )
        shutil.move(src=sim_dir / f"{mode}.tar.gz", dst=mode_compressed)
        for xmm_mode_dir in xmm_filter_dir.rglob(f"{mode}{os.sep}"):
            shutil.rmtree(xmm_mode_dir)
        shutil.rmtree(mode_dir)


def run(path_to_cfg: Path) -> None:
    starttime = datetime.now()
    with open(path_to_cfg) as f:
        cfg: dict[str, dict] = json.load(f)

    global env_cfg, sim_cfg

    env_cfg = EnvironmentCfg(**cfg.pop("environment"))
    sim_cfg = SimulationCfg(
        **cfg.pop("simulation"),
        simput_dir=env_cfg.working_dir / "simput",
        out_dir=env_cfg.working_dir / "xmm_sim_dataset",
    )

    del cfg

    configure_logger(
        log_dir=env_cfg.log_dir,
        log_name=f"03_xmm_simulation_{'_'.join(sim_cfg.instruments)}.log",
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
        for instrument_name in sim_cfg.instruments:
            instrument_dir = xml_dir / instrument_name
            for res_mult in sim_cfg.res_mults:
                res_mult_dir = instrument_dir / sim_cfg.filter / f"{res_mult}x"
                res_mult_dir.mkdir(exist_ok=True, parents=True)

        if env_cfg.working_dir != env_cfg.output_dir:
            simput_compressed = env_cfg.output_dir / "simput.tar.gz"
            if simput_compressed.exists():
                logger.info(f"Found compressed SIMPUT files in {simput_compressed.resolve()}. Decompressing...")
                decompress_targz(
                    in_file_path=simput_compressed,
                    out_file_dir=sim_cfg.simput_dir,
                )
                logger.success("DONE\tDecompressing SIMPUT files.")

        logger.info("START\tCreating all PSF files.")
        to_run = partial(create_psf_file, xml_dir=xml_dir)
        kwds = (
            {"instrument_name": instrument_name, "res_mult": res_mult}
            for instrument_name in sim_cfg.instruments
            for res_mult in sim_cfg.res_mults
        )
        _, duration = mp_run(to_run, kwds, sim_cfg.num_processes, env_cfg.debug)
        logger.success(f"DONE\tPSF files have been created. Duration: {duration}")

        logger.info("START\tCreating all vignetting files.")
        to_run = partial(create_vinget_file, xml_dir=xml_dir)
        kwds = ({"instrument_name": instrument_name} for instrument_name in sim_cfg.instruments)
        _, duration = mp_run(to_run, kwds, sim_cfg.num_processes, env_cfg.debug)
        logger.success(f"DONE\tVignetting files have been created. Duration: {duration}")

        logger.info("START\tCreating all XML files.")
        to_run = partial(
            create_xml_files,
            xml_dir=xml_dir,
            xmm_filter=sim_cfg.filter,
            sim_separate_ccds=sim_cfg.sim_separate_ccds,
            wait_time=sim_cfg.wait_time,
        )
        kwds = (
            {
                "instrument_name": instrument_name,
                "res_mult": res_mult,
            }
            for res_mult in sim_cfg.res_mults
            for instrument_name in sim_cfg.instruments
        )
        _, duration = mp_run(to_run, kwds, sim_cfg.num_processes, env_cfg.debug)
        logger.success(f"DONE\tXML files have been created. Duration: {duration}")

        for instrument_name in sim_cfg.instruments:
            xmm_filter_dir = sim_cfg.out_dir / instrument_name / sim_cfg.filter
            xmm_filter_dir.mkdir(exist_ok=True, parents=True)
            if sim_cfg.modes.img != 0:
                _simulate_mode(
                    instrument_name=instrument_name,
                    mode="img",
                    amount=sim_cfg.modes.img,
                    sim_dir=sim_dir,
                    xmm_filter_dir=xmm_filter_dir,
                    xml_dir=xml_dir,
                )

            if sim_cfg.modes.agn != 0:
                _simulate_mode(
                    instrument_name=instrument_name,
                    mode="agn",
                    amount=sim_cfg.modes.agn,
                    sim_dir=sim_dir,
                    xmm_filter_dir=xmm_filter_dir,
                    xml_dir=xml_dir,
                )

            if sim_cfg.modes.bkg != 0:
                _simulate_mode(
                    instrument_name=instrument_name,
                    mode="bkg",
                    amount=sim_cfg.modes.bkg,
                    sim_dir=sim_dir,
                    xmm_filter_dir=xmm_filter_dir,
                    xml_dir=xml_dir,
                )

    endtime = datetime.now()
    logger.info(f"Duration: {endtime - starttime}")


if __name__ == "__main__":
    parser = ArgumentParser(prog="", description="")
    parser.add_argument("-p", "--config_path", type=Path, required=True, help="Path to config file.")

    args = parser.parse_args()
    run(path_to_cfg=args.config_path)
