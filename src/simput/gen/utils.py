from pathlib import Path
from typing import Union, Tuple, Iterable
from src.utils.external_run import run_headas_command
import numpy as np
from astropy.io import fits


def ones_like_xmm(
        resolution: Union[int, Tuple[int, int]],
        cdelt: float,
        crpix1: int,
        crpix2: int,
        run_dir: Path,
        filename: str
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
        "CDELT1": cdelt,
        "CDELT2": cdelt,
        "comment": "This fits image has all pixel value as 1 and has a similar resolution as xmm"
    }

    header = fits.Header(header)
    hdu = fits.PrimaryHDU(data=np.ones(resolution), header=header)

    out_file = run_dir / filename
    hdu.writeto(out_file, overwrite=True)

    return out_file


def generate_ascii_spectrum(
        run_dir: Path,
        energies,
        rates: Union[float, Iterable],
        verbose: bool
) -> Path:
    if isinstance(rates, float):
        content = [f"{energy} {rates}" for energy in energies]
    elif isinstance(rates, Iterable):
        content = [f"{energy} {rate}" for energy, rate in zip(energies, rates)]
    else:
        raise AttributeError
    content = "\n".join(content)

    out_file = run_dir / "background_spectrum.txt"
    with open(out_file, "w") as f:
        f.write(content)

    if verbose:
        print("energy (keV) | rate (photon/s/cm**2/keV)")
        print(content)
        print(f"Ascii spectrum generated and saved to: {out_file}")

    return out_file


def generate_simput(
        run_dir: Path,
        filename: str,
        emin,
        emax,
        image_file: Path,
        ascii_spectrum_file: Path
) -> Path:
    outfile_path = run_dir / filename

    command = f"simputfile Simput={outfile_path.resolve()} RA=0.0 DEC=0.0 Emin={emin} Emax={emax} " \
              f"history=True clobber=True ImageFile={image_file.resolve()} " \
              f"ASCIIFile={ascii_spectrum_file.resolve()}"

    run_headas_command(command)

    return outfile_path