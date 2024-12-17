import os
import pathlib
import shutil
import tomllib
from argparse import ArgumentParser
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

from dotenv import load_dotenv
from loguru import logger
from tqdm import tqdm

from src.config import DownloadCfg, EnergyCfg, EnvironmentCfg, MultiprocessingCfg
from src.illustris_tng.fits import cutout_to_xray_fits
from src.illustris_tng.web_api import (
    get_available_simulations,
    get_cutouts,
    get_subhalos,
)
from src.tools.files import compress_targz, decompress_targz
from src.tools.run_utils import configure_logger
from src.tools.tools import download_file


def download_data(
    download_cfg: DownloadCfg,
    env_cfg: EnvironmentCfg,
    energies: EnergyCfg,
    mp_cfg: MultiprocessingCfg,
    api_key: str,
) -> None:
    decompress_fs: dict[str, Future] = {}
    compress_fs: dict[Future, str] = {}

    with ProcessPoolExecutor(max_workers=mp_cfg.num_cores) as executor:
        if download_cfg.cutouts_compressed.exists():
            logger.info(f"Found compressed cutouts in {download_cfg.cutouts_compressed}")
            decompress_fs["cutouts"] = executor.submit(
                decompress_targz,
                in_file_path=download_cfg.cutouts_compressed,
                out_file_dir=download_cfg.cutouts_path.parent,
            )
        if download_cfg.fits_compressed.exists():
            logger.info(f"Found compressed FITS in {download_cfg.fits_compressed}")
            decompress_fs["fits"] = executor.submit(
                decompress_targz,
                in_file_path=download_cfg.fits_compressed,
                out_file_dir=download_cfg.fits_path.parent,
            )

        logger.info("START\tGetting simulations.")
        simulations: list[str] = [
            simulation
            for simulation in get_available_simulations(api_key=api_key)
            if simulation[0] in download_cfg.simulations
        ]
        subhalos_fs: list[Future] = []
        for simulation_name, simulation_url in simulations:
            for snapshot in download_cfg.snapshots:
                logger.info(f"START\tGetting subhalo for simulation {simulation_name} and snapshot {snapshot}.")
                subhalos_fs.append(
                    executor.submit(
                        get_subhalos,
                        api_key=api_key,
                        simulation_url=simulation_url,
                        snapshot_num=snapshot,
                        params={
                            "limit": download_cfg.top_n,
                            "primary_flag": 1,
                            "order_by": "-mass_gas",
                        },
                    )
                )

        if not simulations:
            raise ValueError("No simulations found! Please check your config file.")

        logger.success(
            f"DONE\tFound {len(simulations)} available simulations. "
            f"({', '.join([simulation[0] for simulation in simulations])})"
        )
        del simulations  # Cleanup

        if "cutouts" in decompress_fs:
            fs = decompress_fs.pop("cutouts")
            if fs.running():
                logger.info("Waiting for the decompression of the cutouts to finish.")
                fs.result()

        cutouts_fs: list[Future] = []
        with tqdm(total=len(subhalos_fs), desc="Getting subhalos", leave=False) as pbar:
            for future in as_completed(subhalos_fs):
                for subhalo_url in future.result():
                    logger.success(f"DONE\tGot subhalo {subhalo_url} -> Start getting cutouts.")
                    cutouts_fs.append(
                        executor.submit(
                            get_cutouts,
                            subhalo_url=subhalo_url,
                            api_key=api_key,
                            cutout_datafolder=download_cfg.cutouts_path,
                            fail_on_error=env_cfg.fail_on_error,
                        )
                    )
                pbar.update()
            elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])
        del subhalos_fs  # Cleanup

        if len(cutouts_fs) == 0:
            raise ValueError("No subhalos found! Please check your config file.")

        logger.success(f"DONE\tGot {len(cutouts_fs)} subhalos. Duration: {elapsed_time}")

        if "fits" in decompress_fs:
            fs = decompress_fs.pop("fits")
            if fs.running():
                logger.info("Waiting for the decompression of the FITS to finish.")
                fs.result()
        del decompress_fs  # Cleanup

        cloudy_emissivity = env_cfg.working_dir / "cloudy_emissivity_v2.h5"
        if not cloudy_emissivity.exists():
            download_file("http://yt-project.org/data/cloudy_emissivity_v2.h5", cloudy_emissivity)
        fits_fs: dict[Future, Path] = {}
        with tqdm(total=len(cutouts_fs), desc="Gettings cutouts", leave=False) as pbar:
            for future in as_completed(cutouts_fs):
                cutout = future.result()
                path: Path = cutout["file"]
                logger.success(f"Got cutout {path} -> Start generation of X-ray FITS.")
                fs = executor.submit(
                    cutout_to_xray_fits,
                    cutout=path,
                    output_dir=download_cfg.fits_path,
                    sub={"pos_x": float(cutout["x"]), "pos_y": float(cutout["y"]), "pos_z": float(cutout["z"])},
                    mode_dict=download_cfg.modes,
                    cloudy_emissivity_root=cloudy_emissivity.parent,
                    widths=download_cfg.simulations[path.parts[-3]],
                    resolutions=download_cfg.resolutions,
                    energies=energies,
                    redshift=download_cfg.snapshots[int(path.parts[-2])],
                    environment=env_cfg,
                )
                fits_fs[fs] = path
                pbar.update()
            elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])

        logger.success(f"DONE\tGetting cutouts. Duration: {elapsed_time}")
        del cutouts_fs  # Cleanup
        if env_cfg.tar_and_compress:
            logger.info(
                f"Retrieved cutouts will be compressed and moved to {download_cfg.cutouts_compressed.resolve()}. "
                f"Existing file will be overwritten."
            )
            download_cfg.cutouts_compressed.unlink(missing_ok=True)
            fs = executor.submit(
                compress_targz,
                in_path=download_cfg.cutouts_path,
                out_file_path=download_cfg.cutouts_compressed,
            )
            compress_fs[fs] = "cutouts"

        with tqdm(total=len(fits_fs), desc="Creating FITS from cutouts", leave=False) as pbar:
            for future in as_completed(fits_fs):
                cutout: Path = fits_fs[future]
                logger.success(f"Converted {cutout} to X-ray FITS.")
                pbar.update()
            elapsed_time = pbar.format_interval(pbar.format_dict["elapsed"])
        logger.success(f"DONE\tCreated FITS from cutouts. Duration: {elapsed_time}")
        del fits_fs  # Cleanup

        if env_cfg.tar_and_compress:
            logger.info(
                f"Generated FITS will be compressed and moved to {download_cfg.fits_compressed.resolve()}. "
                f"Existing file will be overwritten."
            )
            download_cfg.fits_compressed.unlink(missing_ok=True)
            fs = executor.submit(
                compress_targz,
                in_path=download_cfg.fits_path,
                out_file_path=download_cfg.fits_compressed,
            )
            compress_fs[fs] = "fits"

        if env_cfg.tar_and_compress:
            with tqdm(total=len(compress_fs), desc="Compressing created files", leave=False) as pbar:
                for future in as_completed(compress_fs):
                    files: str = compress_fs[future]
                    if files == "cutouts":
                        logger.success(f"Compressed cutouts --> Start deleting {download_cfg.cutouts_path}")
                        executor.submit(shutil.rmtree, download_cfg.cutouts_path)
                    else:
                        logger.success(f"Compressed FITS --> Start deleting {download_cfg.fits_path}")
                        executor.submit(shutil.rmtree, download_cfg.fits_path)
                    pbar.update()
            del compress_fs  # Cleanup
            logger.info(
                f"Waiting for deletion of {download_cfg.cutouts_path} and {download_cfg.fits_path} to finish..."
            )


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
        help="TNG API key. If you don't have one, create an account at https://www.tng-project.org/data/",
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
    environment = EnvironmentCfg(**cfg.pop("environment"))

    download = DownloadCfg(
        **cfg.pop("download"),
        cutouts_path=environment.working_dir / "cutouts",
        cutouts_compressed=environment.output_dir / "cutouts.tar.gz",
        fits_path=environment.working_dir / "fits",
        fits_compressed=environment.output_dir / "fits.tar.gz",
    )
    energy = EnergyCfg(**cfg.pop("energy"))
    multiprocessing = MultiprocessingCfg(**cfg.pop("multiprocessing"))

    del cfg

    logger.remove()

    configure_logger(
        log_dir=environment.log_dir,
        log_name="01_download_files_{time}.log",
        enqueue=True,
        debug=environment.debug,
        verbose=environment.verbose,
    )

    starttime = datetime.now()
    download_data(
        download_cfg=download,
        env_cfg=environment,
        energies=energy,
        mp_cfg=multiprocessing,
        api_key=args.api_key,
    )
    endtime = datetime.now()
    logger.success(f"DONE\tDuration: {endtime - starttime}")
