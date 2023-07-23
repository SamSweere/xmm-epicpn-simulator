# This code generates the background noise (sky + instrument + particle) based on a background spectrum
from pathlib import Path

import numpy as np
from astropy.io import fits

from src.simput.constants import CDELT, RESOLUTION, PIXEL_SIZE
from src.simput.gen.utils import ones_like_xmm, generate_ascii_spectrum, generate_simput


def get_ascii_spectrum(run_dir: Path, spectrum_file, verbose) -> Path:
    # Open the background spectrum file (sky + instrument + particle)
    with fits.open(spectrum_file) as hdu:
        bin_factor = hdu['SPECTRUM'].header['SPECDELT']
        channels = hdu['SPECTRUM'].data['CHANNEL']
        energies = channels * bin_factor / 1000

        # Calculate the rate based on the counts
        counts = hdu['SPECTRUM'].data['COUNTS'].astype(np.float32)
        rates = counts / float(hdu['SPECTRUM'].header['EXPOSURE'])

    surface = (PIXEL_SIZE * RESOLUTION) ** 2  # cm**2
    cgi_rates = rates / surface  # photon/s/cm**2/keV

    return generate_ascii_spectrum(
        run_dir,
        energies,
        cgi_rates,
        verbose
    )


def generate_background(
        run_dir: Path,
        spectrum_file,
        emin,
        emax,
        verbose
) -> Path:
    image_file = ones_like_xmm(resolution=RESOLUTION, cdelt=CDELT, crpix1=int(RESOLUTION / 2) - 44,
                               crpix2=219, run_dir=run_dir, filename="const_background.fits")

    ascii_spectrum_file = get_ascii_spectrum(run_dir, spectrum_file, verbose)

    outfile_path = generate_simput(
        run_dir=run_dir,
        filename=f"background.fits",
        emin=emin,
        emax=emax,
        image_file=image_file,
        ascii_spectrum_file=ascii_spectrum_file
    )

    if verbose:
        print("Background generation complete")

    return outfile_path
