from datetime import datetime
from pathlib import Path
from typing import List, Tuple
from uuid import uuid4

import numpy as np
from astropy.io import fits

from src.simput.utils import get_spectrumfile
from src.xmm_utils.external_run import run_headas_command


def _prepare_fits_image(
        run_dir: Path,
        img_path: Path,
        zoom: float = 3,
        sigma_b: float = 10,
        offset_x: float = 0,
        offset_y: float = 0,
        fov: float = 0.8
) -> Tuple[Path, float]:
    # fov is in Degrees
    # offset in range -1, 1. Percentage of possible offset, this scales with the zoom level

    with fits.open(img_path) as hdu_in:
        res = hdu_in['XRAY_PHOTON_INTENSITY_0.5_2.0_KEV'].header['NAXIS1']
        data = hdu_in['XRAY_PHOTON_INTENSITY_0.5_2.0_KEV'].data

    max_offset = (res / 3.0) * (1.0 - 1.0 / zoom)  # The maximum pixel offset
    crpix1 = (res / 2.0) + offset_x * max_offset  # Defines the center pixel x
    crpix2 = (res / 2.0) + offset_y * max_offset  # Defines the center pixel y

    # cdelt give the pixel sizes in degrees
    cdelt = (fov / res) * zoom
    cdelt1 = cdelt
    cdelt2 = cdelt

    header = {
        "CUNIT1": "deg",
        "CUNIT2": "deg",
        "CDELT1": -1.0 * cdelt1,
        "CDELT2": cdelt2,
        "CRPIX1": crpix1,
        "CRPIX2": crpix2,
        "CRVAL1": 0.0,
        "CRVAL2": 0.0,
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
        "MTYPE1": "EQPOS",
        "MFORM1": "RA,DEC"
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

    # target_mean = 1000
    # current_mean = np.mean(center_cutout)

    # Making the center a certain brightness. The $flux_{\mu_B}$ and $flux_{\sigma_B}$ where calculated based on a
    # constant distribution. Our sources will not be constant. We want a certain area percentage $x_p$ of the image to
    # have the flux. We therefore have to scale the flux. Note that flux is distributed based on the input image.
    # $P_{counts} = \frac{\sum inner}{\sum all}$
    # $P_{area} = x_{p}^2$
    # $scaling = \frac{P_{area}}/{P_{counts}}$

    sum_inner = np.sum(center_cutout)
    sum_all = np.sum(data)
    p_counts = sum_inner / sum_all
    p_area = box_size_perc ** 2
    scaling = p_area / p_counts
    # We ceil limit the scaling to max 1. Such that if in the unexpected event that the source is not centered
    # at the center the flux will not reach an extreme level and make the SIXTE simulation take an extreme amount of
    # time.
    scaling = min(1.0, scaling)

    # print("Scaling:", scaling)

    # We also increase the flux with the zoom, since in this case less of the whole will be visible and pixels will
    # cover a bigger part of the fov and are therefore dimmer.

    # Calculated based on 50ks exposure of background and constant source with 1e-11 and 2e-11 fluxes
    flux_mu_b = 6.105610561056106e-12  # Flux needed to reach one background
    flux_sigma_b = 7.878743811881188e-12  # Flux needed to reach one sigma
    flux = (flux_mu_b + sigma_b * flux_sigma_b) * zoom * scaling

    # flux = 2e-11
    #
    # baseline_mean = 9.085126612539098e-06
    #
    # norm_scale = baseline_mean/current_mean
    #
    #
    # background_norm = 0.000017
    # sigma_norm = 0.0000385
    # norm = (background_norm + sigma_b * sigma_norm)*zoom*norm_scale

    # Scale the data to uint16 such that is can be compressed better
    max_val = np.iinfo(np.uint16).max
    data = (data / np.max(data)) * max_val
    data = data.astype(np.uint16)

    hdu = fits.PrimaryHDU(data, header=header)

    tmp_output_file = run_dir / f"{uuid4().int}.fits"
    hdu.writeto(tmp_output_file, overwrite=True)

    return tmp_output_file, flux


def simput_image(
        run_dir: Path,
        img_settings: dict,
        keep_files: bool = False,
        verbose: bool = True
) -> List[Path]:
    img_path_in: Path = img_settings['img_path']
    zooms = img_settings['zoom']
    sigmas_b = img_settings['sigma_b']
    offsets_x = img_settings['offset_x']
    offsets_y = img_settings['offset_y']

    emin = 0.5
    emax = 2.0

    # We move the image file to the whole simput
    ra = 0.0
    dec = 0.0

    # Get the spectrum file
    spectrum_file = get_spectrumfile(run_dir=run_dir, verbose=verbose)

    output_files = []

    for zoom, sigma_b, offset_x, offset_y in zip(zooms, sigmas_b, offsets_x, offsets_y):
        img_path, flux = _prepare_fits_image(run_dir, img_path_in, zoom=zoom, sigma_b=sigma_b, offset_x=offset_x,
                                             offset_y=offset_y)

        name = f"{img_path_in.stem}_p0_{emin}ev_p1_{emax}ev_sb_{sigma_b}_zoom_{zoom}_offx_{offset_x}_offy_{offset_y}"
        name = name.replace(".", "_")
        output_file_name = f"{name}.simput"

        output_file = run_dir / output_file_name

        simput_command = (f"simputfile Simput={output_file.resolve()} RA={ra} Dec={dec} "
                          f"XSPECFile={spectrum_file.resolve()} Emin={emin} Emax={emax} "
                          f"ImageFile={img_path.resolve()} srcFlux={flux}")
        run_headas_command(simput_command, verbose=verbose)

        # TODO: rotate simput by a random factor

        # Add specifics to the simput file
        with fits.open(output_file.resolve(), mode="update") as hdu:
            header = hdu['PRIMARY'].header
            header['INPUT'] = (img_path_in.name, "The image file used as input")
            header['ZOOM'] = (zoom, "The amount the image is enlarged")
            header['SIMGMA_B'] = (sigma_b, "Brightness based on the std of 50ks background.")
            header['FLUX'] = (flux, "The flux of the whole image.")
            header['OFFSET_X'] = (offset_x, "Percentage offset of x")
            header['OFFSET_Y'] = (offset_y, "Percentage offset of y")
            header['P0'] = (emin, "Emin")
            header['P1'] = (emax, "Emax")

            header['COMMENT'] = "The image is used as a distribution map for this flux."
            header['COMMENT'] = "All the calibration is done on 50ks."
            header['COMMENT'] = f"Created by Sam Sweere (samsweere@gmail.com) for ESAC at " \
                                f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}"

        output_files.append(output_file.resolve())

    if not keep_files:
        img_path_in.unlink()

    return output_files
