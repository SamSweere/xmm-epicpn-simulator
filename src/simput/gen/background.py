from pathlib import Path

import numpy as np
from astropy.io import fits

from src.simput.constants import CDELT, RESOLUTION, PIXEL_SIZE
from src.simput.gen.utils import ones_like_xmm, generate_ascii_spectrum, generate_simput


def _get_ascii_spectrum(
        run_dir: Path,
        spectrum_file: Path,
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

    surface = (PIXEL_SIZE * RESOLUTION) ** 2  # cm**2
    cgi_rates = rates / surface  # photon/s/cm**2/keV

    ascii_spectrum = generate_ascii_spectrum(run_dir,
                                             energies,
                                             cgi_rates,
                                             verbose)

    return ascii_spectrum


def background(
        run_dir: Path,
        spectrum_file: Path,
        emin: float,
        emax: float,
        suffix=None,
        verbose: bool = True
) -> Path:
    suffix = "" if suffix is None else f"_{suffix}"
    image_file = ones_like_xmm(resolution=RESOLUTION,
                               cdelt=CDELT,
                               crpix1=int(RESOLUTION / 2) - 44,
                               crpix2=219,
                               run_dir=run_dir,
                               filename=f"const_background{suffix}.fits")

    ascii_spectrum_file = _get_ascii_spectrum(run_dir, spectrum_file, verbose)

    outfile_path = generate_simput(run_dir=run_dir,
                                   filename=f"background{suffix}.simput",
                                   emin=emin,
                                   emax=emax,
                                   image_file=image_file,
                                   ascii_spectrum_file=ascii_spectrum_file)

    if verbose:
        print(f"Background generation complete. Saved to {outfile_path.resolve()}")

    return outfile_path
