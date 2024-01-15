import json
import shutil
from argparse import ArgumentParser
from datetime import datetime, timedelta
from functools import partial
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict, List, Optional, Tuple

import numpy as np
from loguru import logger

from src.config import EnergySettings, EnvironmentCfg, SimputCfg
from src.simput.agn import get_fluxes
from src.simput.gen import simput_generate
from src.simput.utils import get_spectrumfile
from src.xmm_utils.file_utils import compress_targz, decompress_targz
from src.xmm_utils.run_utils import configure_logger
from xmm_utils.multiprocessing import mp_run

logger.remove()


def run(
    path_to_cfg: Path, agn_counts_file: Optional[Path], spectrum_dir: Optional[Path]
) -> None:
    starttime = datetime.now()
    with open(path_to_cfg, "r") as f:
        cfg: Dict[str, dict] = json.load(f)
    env_cfg = EnvironmentCfg(**cfg.pop("environment"))
    simput_cfg = SimputCfg(
        **cfg.pop("simput"),
        simput_dir=env_cfg.working_dir / "simput",
        fits_dir=env_cfg.working_dir / "fits",
        fits_compressed=env_cfg.output_dir / "fits.tar.gz",
    )
    energies = EnergySettings(**cfg.pop("energy"))

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

    with TemporaryDirectory(prefix="simput_") as tmp_dir:
        tmp_dir = Path(tmp_dir)

        if simput_cfg.modes.img != 0:
            logger.info(f"START\tGenerating SIMPUT for mode 'img'...")
            img_path = simput_cfg.simput_dir / "img"
            img_path.mkdir(parents=True, exist_ok=True)

            if simput_cfg.fits_compressed.exists():
                logger.info(
                    f"Found compressed FITS files in {simput_cfg.fits_compressed.resolve()}. Decompressing..."
                )
                decompress_targz(
                    in_file_path=simput_cfg.fits_compressed,
                    out_file_dir=simput_cfg.fits_dir,
                )

            to_create: List[Tuple[Path, int]] = []
            rng = np.random.default_rng()

            amount_img = simput_cfg.modes.img
            for in_file in simput_cfg.fits_dir.rglob("*.fits"):
                if amount_img > 0:
                    amount_img = amount_img - 1

                tng_set, snapshot_num = in_file.parts[-2], in_file.parts[-1]
                # Check how many files have already been generated and how many are left to generate
                sample_num = simput_cfg.num_img_sample
                for _ in (img_path / tng_set / snapshot_num).glob(f"{in_file.stem}*"):
                    sample_num = sample_num - 1
                    if sample_num == 0:
                        break

                # None are left to be generated => Skip
                if sample_num == 0:
                    logger.debug(f"Won't generate any images for {in_file.name}.")
                    if env_cfg.consume_data:
                        in_file.unlink()
                    continue

                logger.info(f"Will generate {sample_num} images for {in_file.name}.")
                to_create.append((in_file, sample_num))

                if amount_img == 0:
                    break

            kwds = (
                {
                    "img_settings": {
                        "img_path": in_file.resolve(),
                        "consume_data": env_cfg.consume_data,
                        "zoom": np.round(
                            rng.uniform(
                                low=simput_cfg.zoom_range[0],
                                high=simput_cfg.zoom_range[1],
                                size=num,
                            ),
                            2,
                        ),
                        "sigma_b": np.round(
                            rng.uniform(
                                low=simput_cfg.sigma_b_range[0],
                                high=simput_cfg.sigma_b_range[1],
                                size=num,
                            ),
                            2,
                        ),
                        "offset_x": np.round(
                            rng.normal(
                                loc=-simput_cfg.offset_std,
                                scale=simput_cfg.offset_std,
                                size=num,
                            ),
                            2,
                        ),
                        "offset_y": np.round(
                            rng.normal(
                                loc=-simput_cfg.offset_std,
                                scale=simput_cfg.offset_std,
                                size=num,
                            ),
                            2,
                        ),
                    }
                }
                for in_file, num in to_create
            )

            to_run = partial(
                simput_generate,
                emin=energies.emin,
                emax=energies.emax,
                mode="img",
                tmp_dir=tmp_dir,
                output_dir=img_path,
            )
            _, duration = mp_run(to_run, kwds, simput_cfg.num_processes, env_cfg.debug)
            logger.info(f"DONE\tGenerating SIMPUT for mode 'img'. Duration: {duration}")

            if env_cfg.working_dir != env_cfg.output_dir:
                img_compressed = env_cfg.output_dir / "simput" / "img.tar.gz"
                img_compressed.parent.mkdir(parents=True, exist_ok=True)
                logger.info(
                    f"Compressing IMG SIMPUT files and moving them to {img_compressed.resolve()}."
                    f"Existing file will be overwritten."
                )
                if img_compressed.exists():
                    img_compressed.unlink()
                compress_targz(in_path=img_path, out_file_path=img_compressed)
                shutil.rmtree(img_path)

        if simput_cfg.modes.bkg:
            from src.xmm.utils import get_fov

            if spectrum_dir is None:
                raise FileNotFoundError(f"{spectrum_dir} does not exist!")

            if not spectrum_dir.is_dir():
                raise NotADirectoryError(f"{spectrum_dir} is not a directory!")

            logger.info(f"START\tGenerating SIMPUT for mode 'bkg'...")
            bkg_path = simput_cfg.simput_dir / "bkg"
            img_settings = []
            for instrument_name in simput_cfg.instruments:
                mode_dir = bkg_path / "bkg"
                mode_dir.mkdir(exist_ok=True, parents=True)

                spectrum_name = f"{instrument_name[1]}{instrument_name[-1]}{simput_cfg.filter}ffg_spectrum.fits"
                spectrum_file = spectrum_dir / instrument_name / spectrum_name

                if not spectrum_file.exists():
                    raise FileNotFoundError(f"{spectrum_file} does not exist!")

                if not spectrum_file.is_file():
                    raise FileNotFoundError(f"{spectrum_file} is not a file!")

                fov = get_fov(instrument_name)

                img_settings.append(
                    {
                        "spectrum_file": spectrum_file,
                        "fov": fov,
                        "instrument_name": instrument_name,
                        "output_dir": mode_dir,
                    }
                )

            kwds = (
                {
                    "img_settings": img_setting,
                    "output_dir": img_setting.pop("output_dir"),
                }
                for img_setting in img_settings
            )

            to_run = partial(
                simput_generate,
                emin=energies.emin,
                emax=energies.emax,
                mode="bkg",
                tmp_dir=tmp_dir,
            )
            _, duration = mp_run(to_run, kwds, simput_cfg.num_processes, env_cfg.debug)
            logger.info(f"DONE\tGenerating SIMPUT for mode 'bkg'. Duration: {duration}")

            if env_cfg.working_dir != env_cfg.output_dir:
                bkg_compressed = env_cfg.output_dir / "simput" / "bkg.tar.gz"
                bkg_compressed.parent.mkdir(parents=True, exist_ok=True)
                logger.info(
                    f"Compressing BKG SIMPUT files and moving them to {bkg_compressed.resolve()}."
                    f"Existing file will be overwritten."
                )
                if bkg_compressed.exists():
                    bkg_compressed.unlink()
                compress_targz(in_path=bkg_path, out_file_path=bkg_compressed)
                shutil.rmtree(bkg_path)

        if simput_cfg.modes.agn > 0:
            if agn_counts_file is None:
                raise FileNotFoundError(f"{agn_counts_file} does not exist!")

            if not agn_counts_file.is_file():
                raise FileNotFoundError(f"{agn_counts_file} is not a file!")

            logger.info(f"START\tGenerating SIMPUT for mode 'agn'...")
            agn_path = simput_cfg.simput_dir / "agn"
            agn_path.mkdir(parents=True, exist_ok=True)

            logger.info(f"Will generate {simput_cfg.modes.agn} AGNs.")
            # Get the spectrum file
            spectrum_file = get_spectrumfile(run_dir=tmp_dir, norm=0.001)
            # Get the fluxes from the agn distribution
            fluxes = get_fluxes(agn_counts_file)
            img_settings = {"spectrum_file": spectrum_file, "fluxes": fluxes}
            if env_cfg.debug:
                img_settings["num"] = simput_cfg.modes.agn
                kwds = ({"img_settings": img_settings} for _ in range(1))
            else:
                img_settings["num"] = 1
                kwds = (
                    {"img_settings": img_settings} for _ in range(simput_cfg.modes.agn)
                )

            to_run = partial(
                simput_generate,
                emin=energies.emin,
                emax=energies.emax,
                mode="agn",
                tmp_dir=tmp_dir,
                output_dir=agn_path,
            )
            _, duration = mp_run(to_run, kwds, simput_cfg.num_processes, env_cfg.debug)
            logger.info(f"DONE\tGenerating SIMPUT for mode 'agn'. Duration: {duration}")

            if env_cfg.working_dir != env_cfg.output_dir:
                agn_compressed = env_cfg.output_dir / "simput" / "agn.tar.gz"
                agn_compressed.parent.mkdir(parents=True, exist_ok=True)
                logger.info(
                    f"Compressing AGN SIMPUT files and moving them to {agn_compressed.resolve()}."
                    f"Existing file will be overwritten."
                )
                if agn_compressed.exists():
                    agn_compressed.unlink()
                compress_targz(in_path=agn_path, out_file_path=agn_compressed)
                shutil.rmtree(agn_path)

    endtime = datetime.now()
    logger.info(f"Duration: {endtime - starttime}")


if __name__ == "__main__":
    parser = ArgumentParser(prog="", description="")
    parser.add_argument(
        "-a", "--agn_counts_file", type=Path, help="Path to agn_counts_cgi."
    )
    parser.add_argument(
        "-p", "--config_path", type=Path, required=True, help="Path to config file."
    )
    parser.add_argument(
        "-s", "--spectrum_dir", type=Path, help="Path to spectrum directory."
    )

    args = parser.parse_args()

    run(
        path_to_cfg=args.config_path,
        agn_counts_file=args.agn_counts_file,
        spectrum_dir=args.spectrum_dir,
    )
