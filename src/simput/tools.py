import os
from collections.abc import Iterable
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import numpy as np
from astropy.io import fits
from loguru import logger
from xspec import Model, Xset


def get_spectrumfile(run_dir: Path, norm=0.01) -> Path:
    spectrum_file = run_dir / "spectrum.xcm"
    if not spectrum_file.exists():
        logger.info(f"Spectrum file at {spectrum_file.resolve()} does not exist. Will create a new one.")

        with open(os.devnull, "w") as f, redirect_stdout(f), redirect_stderr(f):
            Model("phabs*power", setPars={1: 0.04, 2: 2.0, 3: norm})
            Xset.save(f"{spectrum_file.resolve()}")

    return spectrum_file


def ones_like_xmm(
    resolution: int | tuple[int, int],
    cdelt1: float,
    cdelt2: float,
    crpix1: float,
    crpix2: float,
    tmp_file: Path,
) -> None:
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

    hdu.writeto(tmp_file, overwrite=True)


def generate_ascii_spectrum(
    ascii_spectrum_file: Path,
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

    content = "\n".join(content)
    content = content.strip()

    with open(ascii_spectrum_file, "w") as f:
        f.write(content)

    logger.info(f"Ascii spectrum generated and saved to: {ascii_spectrum_file}")

    return ascii_spectrum_file
