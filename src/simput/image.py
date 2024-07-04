from datetime import datetime
from pathlib import Path
from uuid import uuid4

import numpy as np
from astropy.io import fits

from src.sixte import commands
from src.tools.files import compress_gzip
from src.xmm.tools import get_fov


def _prepare_fits_image(
    run_dir: Path,
    img_path: Path,
    zoom: float = 3,
    sigma_b: float = 10,
    offset_x: float = 0,
    offset_y: float = 0,
) -> tuple[Path, float]:
    with fits.open(img_path) as hdu_in:
        naxis1 = hdu_in[0].header["NAXIS1"]
        naxis2 = hdu_in[0].header["NAXIS2"]
        data = hdu_in[0].data

    max_offset_x = (naxis1 / 3.0) * (1.0 - 1.0 / zoom)
    max_offset_y = (naxis2 / 3.0) * (1.0 - 1.0 / zoom)
    crpix1 = (naxis1 / 2.0) + offset_x * max_offset_x
    crpix2 = (naxis2 / 2.0) + offset_y * max_offset_y

    # The FOV is the same for EPN, EMOS1, and EMOS2
    fov = get_fov(instrument_name="epn")
    cdelt1 = (fov / naxis1) * zoom
    cdelt2 = (fov / naxis2) * zoom

    header = {
        "CUNIT1": "deg",
        "CUNIT2": "deg",
        "CDELT1": cdelt1,
        "CDELT2": -cdelt2,
        "CRPIX1": crpix1,
        "CRPIX2": crpix2,
        "CRVAL1": 0.0,
        "CRVAL2": 0.0,
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
        "MTYPE1": "EQPOS",
        "MFORM1": "RA,DEC",
    }

    header = fits.Header(header)

    box_size_perc = 0.05
    res = data.shape[0]
    out_pixels = box_size_perc / 2 * res

    center_x = data.shape[0] / 2
    center_y = data.shape[1] / 2
    x_left = int(center_x - out_pixels)
    x_right = int(center_x + out_pixels)
    y_left = int(center_y - out_pixels)
    y_right = int(center_y + out_pixels)

    center_cutout = data[x_left:x_right, y_left:y_right]

    # Making the center a certain brightness. The $flux_{\mu_B}$ and $flux_{\sigma_B}$ where calculated based on a
    # constant distribution. Our sources will not be constant. We want a certain area percentage $x_p$ of the image to
    # have the flux. We therefore have to scale the flux. Note that flux is distributed based on the input image.
    # $P_{counts} = \frac{\sum inner}{\sum all}$
    # $P_{area} = x_{p}^2$
    # $scaling = \frac{P_{area}}/{P_{counts}}$

    sum_inner = np.sum(center_cutout)
    sum_all = np.sum(data)
    p_counts = sum_inner / sum_all
    p_area = box_size_perc**2
    scaling = p_area / p_counts
    # We ceil limit the scaling to max 1. Such that if in the unexpected event that the source is not centered
    # at the center the flux will not reach an extreme level and make the SIXTE simulation take an extreme amount of
    # time.
    # scaling = min(1.0, scaling)

    # We also increase the flux with the zoom, since in this case less of the whole will be visible and pixels will
    # cover a bigger part of the fov and are therefore dimmer.

    # Calculated based on 50ks exposure of background and constant source with 1e-11 and 2e-11 fluxes
    flux_mu_b = 6.105610561056106e-12  # Flux needed to reach one background
    flux_sigma_b = 7.878743811881188e-12  # Flux needed to reach one sigma
    flux = (flux_mu_b + sigma_b * flux_sigma_b) * zoom * scaling

    # Scale the data to uint16 such that is can be compressed better
    max_val = np.iinfo(np.uint16).max
    data = (data / np.max(data)) * max_val
    data = data.astype(np.uint16)

    hdu = fits.PrimaryHDU(data, header=header)

    tmp_output_file = run_dir / f"{uuid4().int}.fits"
    hdu.writeto(tmp_output_file, overwrite=True)

    return tmp_output_file, flux


def simput_image(
    img_path_in: Path,
    emin: float,
    emax: float,
    zooms: np.ndarray,
    sigmas_b: np.ndarray,
    offsets_x: np.ndarray,
    offsets_y: np.ndarray,
    run_dir: Path,
    output_dir: Path,
    xspec_file: Path,
    consume_data: bool,
) -> list[Path]:
    output_files = []
    output_dir.mkdir(parents=True, exist_ok=True)

    for zoom, sigma_b, offset_x, offset_y in zip(zooms, sigmas_b, offsets_x, offsets_y, strict=False):
        img_path, flux = _prepare_fits_image(
            run_dir,
            img_path_in,
            zoom=zoom,
            sigma_b=sigma_b,
            offset_x=offset_x,
            offset_y=offset_y,
        )

        name = f"{img_path_in.stem}_p0_{emin}ev_p1_{emax}ev_sb_{sigma_b}_zoom_{zoom}_offx_{offset_x}_offy_{offset_y}"
        name = name.replace(".", "_")
        output_file_name = f"{name}.simput"

        output_file = run_dir / output_file_name
        compressed_path = output_dir / f"{output_file.name}.gz"

        commands.simputfile(
            simput=output_file,
            ra=0.0,
            dec=0.0,
            src_flux=flux,
            emin=emin,
            emax=emax,
            xspec_file=xspec_file,
            image_file=img_path,
        )

        img_path.unlink()

        # Add specifics to the simput file
        with fits.open(output_file.resolve(), mode="update") as hdu:
            primary_header = hdu["PRIMARY"].header
            primary_header["INPUT"] = (img_path_in.name, "The image file used as input")
            primary_header["ZOOM"] = (zoom, "The amount the image is enlarged")
            primary_header["SIMGMA_B"] = (
                sigma_b,
                "Brightness based on the std of 50ks background.",
            )
            primary_header["FLUX"] = (flux, "The flux of the whole image.")
            primary_header["OFFSET_X"] = (offset_x, "Percentage offset of x")
            primary_header["OFFSET_Y"] = (offset_y, "Percentage offset of y")
            primary_header["P0"] = (emin, "Emin")
            primary_header["P1"] = (emax, "Emax")

            primary_header["COMMENT"] = "The image is used as a distribution map for this flux."
            primary_header["COMMENT"] = "All the calibration is done on 50ks."
            primary_header["COMMENT"] = (
                f"Created by Sam Sweere (samsweere@gmail.com) for ESAC at "
                f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}"
            )

        compress_gzip(output_file, compressed_path, remove_file=True)

        output_files.append(compressed_path.resolve())

    if consume_data:
        img_path_in.unlink()

    return output_files
