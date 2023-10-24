from pathlib import Path
from typing import Literal

import numpy as np
from astropy.io import fits
from loguru import logger

from src.simput.gen.utils import ones_like_xmm, generate_ascii_spectrum, generate_simput
from src.xmm.utils import get_cdelt_for_instrument, get_surface_for_instrument, get_width_height_for_instrument
from src.xmm.utils import get_crpix12_for_instrument


def get_ascii_spectrum(
        run_dir: Path,
        spectrum_file: Path,
        instrument_name: Literal["epn", "emos1", "emos2"],
        verbose: bool = True
) -> Path:
    # Open the background spectrum file (sky + instrument + particle)
    with fits.open(spectrum_file, mode='readonly') as hdu:
        spectrum = hdu['SPECTRUM']
        bin_factor = spectrum.header['SPECDELT']
        channels = spectrum.data['CHANNEL']
        energies = channels * bin_factor / 1000

        # Calculate the rate based on the counts
        counts = spectrum.data['COUNTS'].astype(np.float32)
        rates = counts / float(spectrum.header['EXPOSURE'])

    surface = get_surface_for_instrument(instrument_name=instrument_name, res_mult=1) * 1e-2  # cm**2
    cgi_rates = rates / surface  # photon/s/cm**2/keV

    ascii_spectrum = generate_ascii_spectrum(run_dir,
                                             energies,
                                             cgi_rates,
                                             verbose)

    return ascii_spectrum


def background(
        run_dir: Path,
        spectrum_file: Path,
        instrument_name: Literal["epn", "emos1", "emos2"],
        emin: float,
        emax: float,
        suffix=None,
        verbose: bool = True
) -> Path:
    suffix = "" if suffix is None else f"_{suffix}"

    cdelt = get_cdelt_for_instrument(instrument_name=instrument_name, res_mult=1)
    width, height = get_width_height_for_instrument(instrument_name=instrument_name, res_mult=1)
    crpix1, crpix2 = get_crpix12_for_instrument(instrument_name=instrument_name, res_mult=1)

    image_file = ones_like_xmm(resolution=(width, height),
                               cdelt=cdelt,
                               crpix1=crpix1,
                               crpix2=crpix2,
                               run_dir=run_dir,
                               filename=f"const_background{suffix}.fits")

    ascii_spectrum_file = get_ascii_spectrum(run_dir, spectrum_file, instrument_name, verbose)

    outfile_path = generate_simput(run_dir=run_dir,
                                   filename=f"background{suffix}.simput",
                                   emin=emin,
                                   emax=emax,
                                   image_file=image_file,
                                   ascii_spectrum_file=ascii_spectrum_file)

    if verbose:
        logger.info(f"Background generation complete. Saved to {outfile_path.resolve()}")

    return outfile_path
