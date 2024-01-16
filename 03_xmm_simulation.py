import json
from argparse import ArgumentParser
from datetime import timedelta
from functools import partial
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict

from loguru import logger

from src.sixte.simulator import run_xmm_simulation
from src.xmm.utils import create_psf_file, create_vinget_file, create_xml_files
from src.xmm_utils.multiprocessing import mp_run
from src.xmm_utils.run_utils import configure_logger
from src.xmm_utils.file_utils import compress_targz, decompress_targz
from src.config import EnvironmentCfg, SimulationCfg
from datetime import datetime
import shutil

logger.remove()


def run(path_to_cfg: Path) -> None:
    starttime = datetime.now()
    with open(path_to_cfg, "r") as f:
        cfg: Dict[str, dict] = json.load(f)

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

    with TemporaryDirectory(prefix="xml_") as xml_dir, TemporaryDirectory(
        prefix="sim_"
    ) as sim_dir:
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
                logger.info(
                    f"Found compressed SIMPUT files in {simput_compressed.resolve()}. Decompressing..."
                )
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
        kwds = (
            {"instrument_name": instrument_name}
            for instrument_name in sim_cfg.instruments
        )
        _, duration = mp_run(to_run, kwds, sim_cfg.num_processes, env_cfg.debug)
        logger.success(
            f"DONE\tVignetting files have been created. Duration: {duration}"
        )

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
                logger.info(f"START\tSimulating {instrument_name} for IMG.")
                simputs = (sim_cfg.simput_dir / "img").rglob("*.simput.gz")
                to_run = partial(
                    run_xmm_simulation,
                    instrument_name=instrument_name,
                    xml_dir=xml_dir.resolve(),
                    mode="img",
                    tmp_dir=sim_dir.resolve(),
                    out_dir=xmm_filter_dir.resolve(),
                    exposure=sim_cfg.max_exposure,
                    xmm_filter=sim_cfg.filter,
                    sim_separate_ccds=sim_cfg.sim_separate_ccds,
                    consume_data=env_cfg.consume_data,
                )
                kwds = (
                    {"simput_file": simput.resolve(), "res_mult": res_mult}
                    for res_mult in sim_cfg.res_mults
                    for simput in simputs
                )
                _, duration = mp_run(to_run, kwds, sim_cfg.num_processes, env_cfg.debug)
                logger.success(
                    f"DONE\tSimulating {instrument_name} for IMG. Duration: {duration}"
                )
                if env_cfg.working_dir != env_cfg.output_dir:
                    img_compressed = (
                        env_cfg.output_dir
                        / "xmm_sim_dataset"
                        / instrument_name
                        / sim_cfg.filter
                        / "img.tar.gz"
                    )
                    img_compressed.parent.mkdir(parents=True, exist_ok=True)

                    logger.info(
                        f"Simulated IMG will be compressed and moved to {img_compressed.resolve()}."
                        "Existing file will be overwritten."
                    )

                    compress_targz(
                        in_path=xmm_filter_dir,
                        out_file_path=sim_dir / "img.tar.gz",
                        remove_files=True,
                    )
                    shutil.move(src=sim_dir / "img.tar.gz", dst=img_compressed)
                    for img_dir in xmm_filter_dir.rglob("img"):
                        shutil.rmtree(img_dir)

                if sim_cfg.modes.agn != 0:
                    logger.info(f"START\tSimulating {instrument_name} for AGN.")
                    simputs = (sim_cfg.simput_dir / "agn").rglob("*.simput.gz")
                    to_run = partial(
                        run_xmm_simulation,
                        instrument_name=instrument_name,
                        xml_dir=xml_dir.resolve(),
                        mode="agn",
                        tmp_dir=sim_dir.resolve(),
                        out_dir=xmm_filter_dir.resolve(),
                        exposure=sim_cfg.max_exposure,
                        xmm_filter=sim_cfg.filter,
                        sim_separate_ccds=sim_cfg.sim_separate_ccds,
                        consume_data=env_cfg.consume_data,
                    )
                    kwds = (
                        {"simput_file": simput.resolve(), "res_mult": res_mult}
                        for res_mult in sim_cfg.res_mults
                        for simput in simputs
                    )
                    _, duration = mp_run(
                        to_run, kwds, sim_cfg.num_processes, env_cfg.debug
                    )
                    logger.success(
                        f"DONE\tSimulating {instrument_name} for AGN. Duration: {duration}"
                    )
                    if env_cfg.working_dir != env_cfg.output_dir:
                        agn_compressed = (
                            env_cfg.output_dir
                            / "xmm_sim_dataset"
                            / instrument_name
                            / sim_cfg.filter
                            / "agn.tar.gz"
                        )
                        agn_compressed.parent.mkdir(parents=True, exist_ok=True)

                        logger.info(
                            f"Simulated AGN will be compressed and moved to {agn_compressed.resolve()}."
                            "Existing file will be overwritten."
                        )

                        compress_targz(
                            in_path=xmm_filter_dir,
                            out_file_path=sim_dir / "agn.tar.gz",
                            remove_files=True,
                        )
                        shutil.move(src=sim_dir / "agn.tar.gz", dst=agn_compressed)
                        for agn_dir in xmm_filter_dir.rglob("agn"):
                            shutil.rmtree(agn_dir)

                if sim_cfg.modes.bkg != 0:
                    simput = next(
                        (sim_cfg.simput_dir / "bkg").rglob(
                            f"*{instrument_name}.simput.gz"
                        )
                    )
                    logger.info(f"START\tSimulating {instrument_name} for BKG.")
                    to_run = partial(
                        run_xmm_simulation,
                        instrument_name=instrument_name,
                        xml_dir=xml_dir.resolve(),
                        mode="bkg",
                        tmp_dir=sim_dir.resolve(),
                        out_dir=xmm_filter_dir.resolve(),
                        exposure=sim_cfg.max_exposure,
                        xmm_filter=sim_cfg.filter,
                        sim_separate_ccds=sim_cfg.sim_separate_ccds,
                        consume_data=env_cfg.consume_data,
                    )
                    kwds = (
                        {"simput_file": simput.resolve(), "res_mult": res_mult}
                        for res_mult in sim_cfg.res_mults
                        for _ in range(sim_cfg.modes.bkg)
                    )
                    _, duration = mp_run(
                        to_run, kwds, sim_cfg.num_processes, env_cfg.debug
                    )
                    logger.success(
                        f"DONE\tSimulating {instrument_name} for BKG. Duration: {duration}"
                    )
                    if env_cfg.working_dir != env_cfg.output_dir:
                        bkg_compressed = (
                            env_cfg.output_dir
                            / "xmm_sim_dataset"
                            / instrument_name
                            / sim_cfg.filter
                            / "bkg.tar.gz"
                        )
                        bkg_compressed.parent.mkdir(parents=True, exist_ok=True)

                        logger.info(
                            f"Simulated BKG will be compressed and moved to {bkg_compressed.resolve()}."
                            "Existing file will be overwritten."
                        )

                        compress_targz(
                            in_path=xmm_filter_dir,
                            out_file_path=sim_dir / "bkg.tar.gz",
                            remove_files=True,
                        )
                        shutil.move(src=sim_dir / "bkg.tar.gz", dst=bkg_compressed)
                        for bkg_dir in xmm_filter_dir.rglob("bkg"):
                            shutil.rmtree(bkg_dir)

    endtime = datetime.now()
    logger.info(f"Duration: {endtime - starttime}")


if __name__ == "__main__":
    parser = ArgumentParser(prog="", description="")
    parser.add_argument(
        "-p", "--config_path", type=Path, required=True, help="Path to config file."
    )

    args = parser.parse_args()
    run(path_to_cfg=args.config_path)
