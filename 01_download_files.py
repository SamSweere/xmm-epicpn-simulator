import os
import pathlib
import shutil
import tomllib
from argparse import ArgumentParser
from datetime import datetime, timedelta
from functools import partial
from pathlib import Path

import requests
from dotenv import load_dotenv
from loguru import logger

from src.config import DownloadCfg, EnergySettings, EnvironmentCfg
from src.illustris_tng.download_data import (
    get_available_simulations,
    get_cutouts,
    get_subhalos,
)
from src.illustris_tng.fits import cutout_to_xray_fits
from src.xmm_utils.file_utils import compress_targz, decompress_targz
from src.xmm_utils.multiprocessing import mp_run
from src.xmm_utils.run_utils import configure_logger

logger.remove()


def run(path_to_cfg: Path, api_key: str):
    starttime = datetime.now()
    with open(path_to_cfg, "rb") as file:
        cfg: dict[str, dict] = tomllib.load(file)
    env_cfg = EnvironmentCfg(**cfg.pop("environment"))

    download_cfg = DownloadCfg(
        **cfg.pop("download"),
        cutouts_path=env_cfg.working_dir / "cutouts",
        cutouts_compressed=env_cfg.output_dir / "cutouts.tar.gz",
        fits_path=env_cfg.working_dir / "fits",
        fits_compressed=env_cfg.output_dir / "fits.tar.gz",
    )
    energies = EnergySettings(**cfg.pop("energy"))

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

    logger.info("START\tGetting simulations.")
    simulations: list[str] = []
    for simulation_name, simulation_url in get_available_simulations(api_key=api_key):
        if simulation_name in download_cfg.simulations:
            simulations.append(simulation_url)

    if not simulations:
        raise ValueError(f"No simulations found! Please check your config file ({path_to_cfg}).")

    logger.info(f"DONE\tFound {len(simulations)} available simulations. ({', '.join(simulations)})")

    logger.info("START\tGetting subhalos.")
    subhalos = []
    for url in simulations:
        for snapshot_num in download_cfg.snapshots:
            # request and inspect most massive top_n subhalos that are central (primary_flag = 1)
            # primary_flag = 1 indicates that this is the central (i.e. most massive, or "primary") subhalo of
            # this FoF halo.
            # see https://www.tng-project.org/data/docs/specifications/#sec2a
            urls = get_subhalos(
                api_key=api_key,
                simulation_url=url,
                snapshot_num=snapshot_num,
                params={
                    "limit": download_cfg.top_n,
                    "primary_flag": 1,
                    "order_by": "-mass_gas",
                },
            )
            for subhalo_url in urls:
                subhalos.append(subhalo_url)

    if not subhalos:
        raise ValueError(f"No subhalos found! Please check your config file ({path_to_cfg}).")

    logger.info(f"DONE\tGot {len(subhalos)} subhalos.")

    logger.info("START\tDownloading/getting already downloaded cutouts.")
    if download_cfg.cutouts_compressed.exists():
        logger.info(f"Found compressed cutouts in {download_cfg.cutouts_compressed.resolve()}. Decompressing...")
        decompress_targz(
            in_file_path=download_cfg.cutouts_compressed,
            out_file_dir=download_cfg.cutouts_path,
        )

    to_run = partial(
        get_cutouts,
        api_key=api_key,
        cutout_datafolder=download_cfg.cutouts_path,
        fail_on_error=env_cfg.fail_on_error,
    )

    kwds = ({"subhalo_url": subhalo} for subhalo in subhalos)
    cutouts, duration = mp_run(to_run, kwds, download_cfg.num_processes, env_cfg.debug)
    logger.info(f"DONE\tDownloading/getting cutouts. Duration: {duration}")

    if env_cfg.working_dir != env_cfg.output_dir:
        logger.info(
            f"Downloaded cutouts will be compressed and moved to {download_cfg.cutouts_compressed.resolve()}."
            f"Existing file will be overwritten."
        )
        if download_cfg.cutouts_compressed.exists():
            download_cfg.cutouts_compressed.unlink()
        compress_targz(
            in_path=download_cfg.cutouts_path,
            out_file_path=download_cfg.cutouts_compressed,
        )

    logger.info("START\tGenerating FITS from cutouts.")
    cloudy_emissivity = env_cfg.working_dir / "cloudy_emissivity_v2.h5"
    if not cloudy_emissivity.exists():
        logger.info(f"Downloading cloudy_emissivity_v2.h5 to {cloudy_emissivity.resolve()}")
        retries = 3
        while retries > 0:
            try:
                with requests.get(
                    "http://yt-project.org/data/cloudy_emissivity_v2.h5",
                    stream=True,
                ) as r:
                    r.raise_for_status()
                    with open(cloudy_emissivity, "wb") as f:
                        for chunk in r.iter_content(chunk_size=int(1e6)):
                            f.write(chunk)
                retries = 0
            except:  # noqa
                retries = retries - 1

        if not cloudy_emissivity.exists():
            raise FileNotFoundError(f"Failed to load cloudy_emissivity_v2.h5 {cloudy_emissivity}!")

    if download_cfg.fits_compressed.exists():
        logger.info(f"Found compressed FITS in {download_cfg.fits_compressed.resolve()}. Decompressing...")
        decompress_targz(
            in_file_path=download_cfg.fits_compressed,
            out_file_dir=download_cfg.fits_path,
        )

    to_run = partial(
        cutout_to_xray_fits,
        output_dir=download_cfg.fits_path,
        mode_dict=download_cfg.modes,
        cloudy_emissivity_root=cloudy_emissivity.parent,
        emin=energies.emin,
        emax=energies.emax,
        resolutions=download_cfg.resolutions,
        overwrite=env_cfg.overwrite,
        fail_on_error=env_cfg.fail_on_error,
        consume_data=env_cfg.consume_data,
    )

    logger.info(f"FITS images will be generated for {len(cutouts)} cutouts.")
    kwds = (
        {
            "cutout": cutout["file"],
            "sub": {
                "pos_x": float(cutout["x"]),
                "pos_y": float(cutout["y"]),
                "pos_z": float(cutout["z"]),
            },
            "widths": download_cfg.simulations[cutout["file"].parts[-3]],
            "redshift": download_cfg.snapshots[int(cutout["file"].parts[-2])],
        }
        for cutout in cutouts
    )
    _, duration = mp_run(to_run, kwds, download_cfg.num_processes, env_cfg.debug)
    logger.info(f"DONE\tGenerating FITS from cutouts. Duration: {duration}")

    if env_cfg.working_dir != env_cfg.output_dir:
        logger.info(
            f"Generated FITS will be compressed and moved to {download_cfg.fits_compressed.resolve()}."
            f"Existing file will be overwritten."
        )
        if download_cfg.fits_compressed.exists():
            download_cfg.fits_compressed.unlink()
        compress_targz(
            in_path=download_cfg.fits_path,
            out_file_path=download_cfg.fits_compressed,
        )

        logger.info(f"Deleting {download_cfg.cutouts_path.resolve()} and {download_cfg.fits_path.resolve()}.")
        shutil.rmtree(download_cfg.cutouts_path)
        shutil.rmtree(download_cfg.fits_path)

    endtime = datetime.now()
    logger.info(f"Duration: {endtime - starttime}")


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
    run(
        path_to_cfg=args.config_path,
        api_key=args.api_key,
    )
