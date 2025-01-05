from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory
from uuid import uuid4

import numpy as np
from astropy.io import fits

import src.heasoft as hsp
from src.config import EnergyCfg, SimputImgCfg
from src.sixte import commands
from src.xmm.tools import get_fov


def _prepare_fits_image(
    tmp_file: Path,
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

    hdu.writeto(tmp_file, overwrite=True)

    return flux


def simput_image(
    img_path_in: Path,
    energies: EnergyCfg,
    amount: int,
    cfg: SimputImgCfg,
    output_dir: Path,
    xspec_file: Path,
    consume_data: bool,
) -> list[Path]:
    output_files = []
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng()

    zooms = np.round(
        rng.uniform(
            cfg.zoom_range[0],
            cfg.zoom_range[1],
            amount,
        ),
        2,
    )
    sigmas_b = np.round(
        rng.uniform(
            cfg.sigma_b_range[0],
            cfg.sigma_b_range[1],
            amount,
        ),
        2,
    )

    offsets_x = np.round(
        rng.normal(
            -cfg.offset_std,
            cfg.offset_std,
            amount,
        ),
        2,
    )
    offsets_y = np.round(
        rng.normal(
            -cfg.offset_std,
            cfg.offset_std,
            amount,
        ),
        2,
    )

    with TemporaryDirectory(prefix="simput_img_") as run_dir:
        run_dir = Path(run_dir)
        for zoom, sigma_b, offset_x, offset_y in zip(zooms, sigmas_b, offsets_x, offsets_y, strict=False):
            out_path = output_dir / f"{img_path_in.stem}_{uuid4().int}.simput.gz"

            with (
                NamedTemporaryFile(mode="r", dir=run_dir, suffix=".fits") as image_file,
                NamedTemporaryFile(mode="r", dir=run_dir, suffix=".simput") as local_out,
                NamedTemporaryFile(mode="w", dir=run_dir) as tmp_file,
            ):
                flux = _prepare_fits_image(
                    tmp_file=image_file.name,
                    img_path=img_path_in,
                    zoom=zoom,
                    sigma_b=sigma_b,
                    offset_x=offset_x,
                    offset_y=offset_y,
                )

                commands.simputfile(
                    simput=local_out.name,
                    ra=0.0,
                    dec=0.0,
                    src_flux=flux,
                    emin=energies.emin,
                    emax=energies.emax,
                    xspec_file=xspec_file,
                    image_file=image_file.name,
                )

                # Add specifics to the simput file
                tmp_file.write(f"INPUT = {img_path_in.name} / The image file used as input\n")
                tmp_file.write(f"ZOOM = {zoom} / The amount the image is enlarged\n")
                tmp_file.write(f"SIGMA_B = {sigma_b} / Brightness based on the std of 50ks background\n")
                tmp_file.write(f"FLUX = {flux} / The flux of the whole image\n")
                tmp_file.write(f"OFFSET_X = {offset_x} / Percentage offset of x\n")
                tmp_file.write(f"OFFSET_Y = {offset_y} / Percentage offset of y\n")
                tmp_file.write(f"P0 = {energies.emin} / Emin\n")
                tmp_file.write(f"P1 = {energies.emax} / Emax\n")
                tmp_file.write("COMMENT = The image is used as a distribution map for this flux.\n")
                tmp_file.write("COMMENT = All the calibration is done on 50ks.\n")
                tmp_file.write("COMMENT = The image is used as a distribution map for this flux.")
                tmp_file.seek(0)

                hsp.fthedit(infile=f"{local_out.name}", keyword=f"@{tmp_file.name}")
                hsp.ftcopy(infile=f"{local_out.name}", outfile=f"{out_path}")
            output_files.append(out_path)

    if consume_data:
        img_path_in.unlink()

    return output_files
