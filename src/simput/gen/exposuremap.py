from pathlib import Path
from typing import Literal

from astropy.io import fits
from loguru import logger

from src.simput.gen.utils import ones_like_xmm, generate_ascii_spectrum, generate_simput
from src.xmm.utils import get_cdelt_for_instrument, get_width_height_for_instrument, get_crpix12_for_instrument


def get_ascii_spectrum(
        run_dir: Path,
        spectrum_file: Path,
        verbose: bool = True
) -> Path:
    # Open the background spectrum file (sky + instrument + particle)
    with fits.open(spectrum_file, mode='readonly') as hdu:
        bin_factor = hdu['SPECTRUM'].header['SPECDELT']
        channels = hdu['SPECTRUM'].data['CHANNEL']

    energies = channels * bin_factor / 1000
    rate = 1.0e-1

    ascii_spectrum = generate_ascii_spectrum(run_dir,
                                             energies,
                                             rate,
                                             verbose)

    return ascii_spectrum


def exposure_map(
        run_dir: Path,
        spectrum_file: Path,
        instrument_name: Literal["epn", "emos1", "emos2"],
        emin: float,
        emax: float,
        suffix=None,
        verbose: bool = True
):
    suffix = "" if suffix is None else f"_{suffix}"

    cdelt = get_cdelt_for_instrument(instrument_name=instrument_name, res_mult=1)
    width, height = get_width_height_for_instrument(instrument_name=instrument_name, res_mult=1)
    crpix1, crpix2 = get_crpix12_for_instrument(instrument_name=instrument_name, res_mult=1)

    image_file = ones_like_xmm(resolution=(width, height),
                               cdelt=cdelt,
                               crpix1=crpix1,
                               crpix2=crpix2,
                               run_dir=run_dir,
                               filename=f"exposure_map{suffix}.fits")

    ascii_spectrum_file = get_ascii_spectrum(run_dir, spectrum_file, verbose)

    outfile_path = generate_simput(run_dir=run_dir,
                                   filename=f"exposure_map{suffix}.simput",
                                   emin=emin,
                                   emax=emax,
                                   image_file=image_file,
                                   ascii_spectrum_file=ascii_spectrum_file)

    if verbose:
        logger.info(f"Exposure map generation complete. Saved to {outfile_path.resolve()}")

    return outfile_path
