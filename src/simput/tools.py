import os
import shutil
from collections.abc import Iterable
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import numpy as np
from astropy.io import fits
from loguru import logger
from xspec import Model, Xset

from src.sixte import commands


def get_spectrumfile(run_dir: Path, norm=0.01) -> Path:
    spectrum_file = run_dir / "spectrum.xcm"
    if not spectrum_file.exists():
        logger.info(f"Spectrum file at {spectrum_file.resolve()} does not exist. Will create a new one.")

        with open(os.devnull, "w") as f, redirect_stdout(f), redirect_stderr(f):
            Model("phabs*power", setPars={1: 0.04, 2: 2.0, 3: norm})
            Xset.save(f"{spectrum_file.resolve()}")

    return spectrum_file


def merge_simputs(simput_files: list[Path], output_file: Path) -> Path:
    # Combine the simput point sources
    if len(simput_files) == 1:
        file = simput_files[0]
        shutil.copy2(file, output_file)
    else:
        commands.simputmerge(infiles=simput_files, outfile=output_file, fetch_extension=True)

    return output_file


def ones_like_xmm(
    resolution: int | tuple[int, int],
    cdelt1: float,
    cdelt2: float,
    crpix1: float,
    crpix2: float,
    run_dir: Path,
    filename: str,
) -> Path:
    if isinstance(resolution, int):
        resolution = (resolution, resolution)

    header = {
        "MTYPE1": "EQPOS",
        "MFORM1": "RA,DEC",
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
        "CRPIX1": crpix1,
        "CRPIX2": crpix2,
        "CRVAL1": 0.0,
        "CRVAL2": 0.0,
        "CUNIT1": "deg",
        "CUNIT2": "deg",
        "CDELT1": cdelt1,
        "CDELT2": cdelt2,
        "comment": "This fits image has all pixel value as 1 and has a similar resolution as xmm",
    }

    header = fits.Header(header)
    hdu = fits.PrimaryHDU(data=np.ones(resolution), header=header)

    out_file = run_dir / filename
    hdu.writeto(out_file, overwrite=True)

    return out_file


def generate_ascii_spectrum(
    run_dir: Path,
    energies: float | Iterable | np.ndarray,
    rates: float | Iterable | np.ndarray,
) -> Path:
    if not isinstance(energies, float | Iterable | np.ndarray):
        raise TypeError(f"'energies' has to be one of (float, Iterable, np.ndarray)! Got: {type(energies)}")
    if not isinstance(rates, float | Iterable | np.ndarray):
        raise TypeError(f"'rates' has to be one of (float, Iterable, np.ndarray)! Got: {type(rates)}")

    if isinstance(energies, float):
        if not isinstance(rates, float):
            raise ValueError("If 'energies' is a float, than 'rates' has to be a float too!")
        content = [f"{energies} {rates}"]
    else:
        rates = rates if isinstance(rates, Iterable) else [rates for _ in energies]
        content = [f"{energy} {rate}" for energy, rate in zip(energies, rates, strict=False)]

    content = f"{os.linesep}".join(content)
    content = content.strip()

    out_file = run_dir / "ascii_spectrum.txt"
    with open(out_file, "w") as f:
        f.write(content)

    logger.info(f"Ascii spectrum generated and saved to: {out_file.resolve()}")

    return out_file
