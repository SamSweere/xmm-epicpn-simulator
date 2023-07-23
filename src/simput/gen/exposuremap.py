# This code generates the background noise (sky + instrument + particle) based on a background spectrum
from pathlib import Path

from astropy.io import fits

from src.simput.constants import CDELT, RESOLUTION
from src.simput.gen.utils import ones_like_xmm, generate_ascii_spectrum, generate_simput


def get_ascii_spectrum(run_dir: Path, spectrum_file, verbose):
    # Open the background spectrum file (sky + instrument + particle)
    with fits.open(spectrum_file) as hdu:
        bin_factor = hdu['SPECTRUM'].header['SPECDELT']
        channels = hdu['SPECTRUM'].data['CHANNEL']

    energies = channels * bin_factor / 1000
    rate = 1.0e-1

    return generate_ascii_spectrum(
        run_dir,
        energies,
        rate,
        verbose
    )


def generate_exposure_map(run_dir: Path, spectrum_file, emin, emax, verbose):
    image_file = ones_like_xmm(resolution=RESOLUTION, cdelt=CDELT, crpix1=186, crpix2=219, run_dir=run_dir,
                               filename="exposure_map.fits")

    ascii_spectrum_file = get_ascii_spectrum(run_dir, spectrum_file, verbose)

    outfile_path = generate_simput(
        run_dir=run_dir,
        filename="exposure_map.simput",
        emin=emin,
        emax=emax,
        image_file=image_file,
        ascii_spectrum_file=ascii_spectrum_file
    )

    if verbose:
        print("Exposure map generation complete")

    return outfile_path
