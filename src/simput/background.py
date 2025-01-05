from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Literal

import numpy as np
from astropy.io import fits
from loguru import logger

import src.heasoft as hsp
from src.config import EnergyCfg
from src.simput.tools import generate_ascii_spectrum, ones_like_xmm
from src.sixte import commands
from src.xmm.tools import get_cdelt, get_crpix12, get_naxis12, get_pixel_size


def get_ascii_spectrum(
    ascii_spectrum_file: Path,
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

    return generate_ascii_spectrum(ascii_spectrum_file, energies, cgi_rates)


def create_background(
    output_dir: Path,
    spectrum_file: Path,
    instrument_name: Literal["epn", "emos1", "emos2"],
    energies: EnergyCfg,
) -> Path:
    suffix = f"_{instrument_name}_{energies.emin}keV_{energies.emax}keV"

    cdelt1, cdelt2 = get_cdelt(instrument_name=instrument_name, res_mult=1)
    naxis1, naxis2 = get_naxis12(instrument_name=instrument_name, res_mult=1)
    crpix1, crpix2 = get_crpix12(instrument_name, 1)

    surface = (get_pixel_size(instrument_name, 1) ** 2) * naxis1 * naxis2 * 1e-2  # cm**2

    outfile = output_dir / f"background{suffix}.simput.gz"

    with (
        NamedTemporaryFile(mode="r", prefix="bkg_", suffix=".simput") as local_out,
        NamedTemporaryFile(mode="r", prefix="bkg_", suffix=".fits") as image_file,
        NamedTemporaryFile(mode="r", prefix="asci_spectrum_bkg", suffix=".txt") as ascii_spectrum_file,
    ):
        ones_like_xmm(
            resolution=(naxis1, naxis2),
            cdelt1=cdelt1,
            cdelt2=cdelt2,
            crpix1=crpix1,
            crpix2=crpix2,
            tmp_file=image_file.name,
        )

        ascii_spectrum_file = get_ascii_spectrum(Path(ascii_spectrum_file.name), spectrum_file, surface)

        commands.simputfile(
            simput=local_out.name,
            emin=energies.emin,
            emax=energies.emax,
            ascii_file=ascii_spectrum_file,
            image_file=image_file.name,
        )

        hsp.ftcopy(infile=local_out.name, outfile=outfile)

    logger.info(f"Background generation complete. Saved to {outfile}")

    return outfile
