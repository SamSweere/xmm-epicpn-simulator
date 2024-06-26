import os
from collections.abc import Iterable
from pathlib import Path

import numpy as np
from astropy.io import fits
from loguru import logger


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
    verbose: bool = True,
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

    if verbose:
        logger.info(f"Ascii spectrum generated and saved to: {out_file.resolve()}")

    return out_file
