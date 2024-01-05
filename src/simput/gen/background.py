from pathlib import Path
from typing import Literal

import numpy as np
from astropy.io import fits
from loguru import logger

from src.simput.gen.utils import ones_like_xmm, generate_ascii_spectrum
from src.sixte import commands
from src.xmm.utils import get_cdelt, get_surface, get_naxis12


def get_ascii_spectrum(
    run_dir: Path,
    spectrum_file: Path,
    instrument_name: Literal["epn", "emos1", "emos2"],
    verbose: bool = True,
) -> Path:
    # Open the background spectrum file (sky + instrument + particle)
    with fits.open(spectrum_file, mode="readonly") as hdu:
        spectrum = hdu["SPECTRUM"]
        bin_factor = spectrum.header["SPECDELT"]
        channels = spectrum.data["CHANNEL"]
        energies = channels * bin_factor / 1000

        # Calculate the rate based on the counts
        counts = spectrum.data["COUNTS"].astype(np.float32)
        rates = counts / float(spectrum.header["EXPOSURE"])

    surface = get_surface(instrument_name=instrument_name, res_mult=1) * 1e-2  # cm**2
    cgi_rates = rates / surface  # photon/s/cm**2/keV

    ascii_spectrum = generate_ascii_spectrum(run_dir, energies, cgi_rates, verbose)

    return ascii_spectrum


def background(
    run_dir: Path,
    spectrum_file: Path,
    instrument_name: Literal["epn", "emos1", "emos2"],
    emin: float,
    emax: float,
    suffix=None,
    verbose: bool = True,
) -> Path:
    suffix = f"_{instrument_name}" if suffix is None else f"_{suffix}"

    cdelt = get_cdelt(instrument_name=instrument_name, res_mult=1)
    naxis1, naxis2 = get_naxis12(instrument_name=instrument_name, res_mult=1)

    crpix1 = round(((naxis1 + 1) / 2.0), 6)
    crpix2 = round(((naxis1 + 1) / 2.0), 6)

    image_file = ones_like_xmm(
        resolution=(naxis1, naxis2),
        cdelt=cdelt,
        crpix1=crpix1,
        crpix2=crpix2,
        run_dir=run_dir,
        filename=f"const_background{suffix}.fits",
    )

    ascii_spectrum_file = get_ascii_spectrum(
        run_dir, spectrum_file, instrument_name, verbose
    )

    outfile_path = run_dir / f"background{suffix}.simput"

    commands.simputfile(
        simput=outfile_path,
        emin=emin,
        emax=emax,
        ascii_file=ascii_spectrum_file,
        image_file=image_file,
    )

    if verbose:
        logger.info(
            f"Background generation complete. Saved to {outfile_path.resolve()}"
        )

    return outfile_path
