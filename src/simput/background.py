from pathlib import Path
from typing import Literal

import numpy as np
from astropy.io import fits
from loguru import logger

from src.simput.tools import generate_ascii_spectrum, ones_like_xmm
from src.sixte import commands
from src.tools.files import compress_gzip
from src.xmm.tools import get_cdelt, get_crpix12, get_naxis12, get_pixel_size


def get_ascii_spectrum(
    run_dir: Path,
    spectrum_file: Path,
    surface: float,
) -> Path:
    # Open the background spectrum file (sky + instrument + particle)
    spectrum, header = fits.getdata(spectrum_file, "SPECTRUM", header=True)
    bin_factor = header["SPECDElT"]
    channels = spectrum["CHANNEL"]
    counts = spectrum["COUNTS"].astype(np.float32)
    exposure = header["EXPOSURE"]

    energies = channels * bin_factor / 1000
    rates = counts / float(exposure)

    cgi_rates = rates / surface  # photon/s/cm**2/keV

    ascii_spectrum = generate_ascii_spectrum(run_dir, energies, cgi_rates)

    return ascii_spectrum


def create_background(
    run_dir: Path,
    output_dir: Path,
    spectrum_file: Path,
    instrument_name: Literal["epn", "emos1", "emos2"],
    emin: float,
    emax: float,
) -> list[Path]:
    suffix = f"_{instrument_name}_{emin}keV_{emax}keV"

    cdelt1, cdelt2 = get_cdelt(instrument_name=instrument_name, res_mult=1)
    naxis1, naxis2 = get_naxis12(instrument_name=instrument_name, res_mult=1)
    crpix1, crpix2 = get_crpix12(instrument_name, 1)

    image_file = ones_like_xmm(
        resolution=(naxis1, naxis2),
        cdelt1=cdelt1,
        cdelt2=cdelt2,
        crpix1=crpix1,
        crpix2=crpix2,
        run_dir=run_dir,
        filename=f"const_background{suffix}.fits",
    )

    surface = (get_pixel_size(instrument_name, 1) ** 2) * naxis1 * naxis2 * 1e-2  # cm**2

    ascii_spectrum_file = get_ascii_spectrum(run_dir, spectrum_file, surface)

    outfile_path = run_dir / f"background{suffix}.simput"
    compressed_path = output_dir / f"{outfile_path.name}.gz"

    commands.simputfile(
        simput=outfile_path,
        emin=emin,
        emax=emax,
        ascii_file=ascii_spectrum_file,
        image_file=image_file,
    )

    compress_gzip(outfile_path, compressed_path, remove_file=True)

    logger.info(f"Background generation complete. Saved to {compressed_path}")

    return [compressed_path]
