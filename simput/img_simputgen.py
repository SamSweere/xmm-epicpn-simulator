from datetime import datetime
import os
import time

from astropy.io import fits

from simput.spectrum import get_spectrumfile
from simput.utils import copy_to_save_folder
import numpy as np

from utils.external_run import run_headas_command


def prepare_fits_image(run_dir, img_path, zoom=3, sigma_b=10, offset_x=0, offset_y=0, fov=0.8):
    # fov is in Degrees
    # offset in range -1, 1. Percentage of possible offset, this scales with the zoom level

    hdu_in = fits.open(img_path)

    res = hdu_in['XRAY_PHOTON_INTENSITY_0.5_2.0_KEV'].header['NAXIS1']

    maxq_offset = (res / 3.0) * (1.0 - 1.0 / zoom)  # The maximum pixel offset
    crpix1 = (res / 2.0) + offset_x * max_offset  # Defines the center pixel x
    crpix2 = (res / 2.0) + offset_y * max_offset  # Defines the center pixel y

    # cdelt give the pixel sizes in degrees
    cdelt = (fov / res) * zoom
    cdelt1 = cdelt
    cdelt2 = cdelt

    header = fits.Header()

    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CDELT1'] = -1.0 * cdelt1
    header['CDELT2'] = cdelt2
    header['CRPIX1'] = crpix1
    header['CRPIX2'] = crpix2
    header['CRVAL1'] = 0.0
    header['CRVAL2'] = 0.0

    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    header['MTYPE1'] = 'EQPOS'
    header['MFORM1'] = 'RA,DEC'

    data = hdu_in['XRAY_PHOTON_INTENSITY_0.5_2.0_KEV'].data

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
    current_mean = np.mean(center_cutout)

    # Making the center a certain brightness. The $flux_{\mu_B}$ and $flux_{\sigma_B}$ where calculated based on a constant distribution. Our sources will not be constant. We want a certain area percentage $x_p$ of the image to have the flux. We therefore have to scale the flux. Note that flux is distributed based on the input image.
    # $P_{counts} = \frac{\sum inner}{\sum all}$
    # $P_{area} = x_{p}^2$
    # $scaling = \frac{P_{area}}/{P_{counts}}$

    sum_inner = np.sum(center_cutout)
    sum_all = np.sum(data)
    p_counts = sum_inner / sum_all
    p_area = box_size_perc ** 2
    scaling = p_area / p_counts
    # We ceil limit the scaling to max 1. Such that if in the unexpected event that the source is not centered at the center the flux will not reach an extreme level and make the SIXTE simulation take an extreme amount of time.
    scaling = min(1.0, scaling)

    print("Scaling:", scaling)

    # We also increase the flux with the zoom, since in this case less of the whole will be visible and pixels will cover a bigger part of the fov and are therefore dimmer.

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

    timestamp = str(time.time()).replace(".", "")[-10:-1]

    tmp_output_file = os.path.join(run_dir, f"tmp_altered_{timestamp}.fits")
    hdu.writeto(tmp_output_file, overwrite=True)

    return tmp_output_file, flux


def create_img_simput(run_dir, img_settings, verbose=True):
    # zoom = 1, sigma_b = 10, offset_x = 0.0, offset_y = 0.0,
    img_path_in = img_settings['img_path']
    zoom = img_settings['zoom']
    sigma_b = img_settings['sigma_b']
    offset_x = img_settings['offset_x']
    offset_y = img_settings['offset_y']

    # open the input fits

    img_name = os.path.split(img_path_in)[-1]
    # sigma_b = 10
    # zoom = 1

    img_path, flux = prepare_fits_image(run_dir, img_path_in, zoom=zoom, sigma_b=sigma_b, offset_x=offset_x,
                                        offset_y=offset_y)

    # simput_files = []
    emin = 0.5
    emax = 2.0

    # We move the image file to the whole simput
    ra = 0.0
    dec = 0.0

    print("flux:", flux)

    # Use the current time as id, such that clashes don't happen
    # id = str(time.time()).replace(".", "")
    fits_name_path = img_name.replace(".fits", "")
    name = f"{fits_name_path}_p0_{emin}ev_p1_{emax}ev_sb_{sigma_b}_zoom_{zoom}_offx_{offset_x}_offy_{offset_y}"
    name = name.replace(".", "_")
    simput_file_name = name + ".simput"
    # timestamp = str(time.time()).replace(".", "")[-10:-1]

    # simput_file_name = f"sim_in_img_{timestamp}.simput"
    output_file_path = os.path.join(run_dir, simput_file_name)

    # Get the spectrum file
    spectrum_file = get_spectrumfile(run_dir=run_dir, verbose=verbose)

    # # Epic pn has a fov of 30 arcmin = 0.5 degrees.
    # pn_fov = 0.5
    #
    # # Randomly position the point source within the fov
    # center_point = (0.0, 0.0)
    #
    # offset = (0.0, 0.0)
    # #TODO: randomize offset
    # # if offset == 'random':
    # #     offset = np.random.uniform(low=-1.0 * pn_fov / 2, high=pn_fov / 2, size=2)
    #
    # location = (center_point[0] + offset[0], center_point[1] + offset[1])
    # ra = location[0]

    # Create a tmp simput name such that the simputfile does not crash
    tmp_simput_path = os.path.join(run_dir, "tmp_img.simput")

    # We need the xspec from headas
    simput_command = f"simputfile RA={ra} Dec={dec} XSPECFile={spectrum_file} Emin={emin} Emax={emax} " \
                     f"Simput={tmp_simput_path} ImageFile={img_path} srcFlux={flux}"
    run_headas_command(simput_command, verbose=verbose)

    # TODO: rotate simput by a random factor

    # Add specifics to the simput file
    hdu = fits.open(tmp_simput_path)
    header = hdu['PRIMARY'].header
    header['INPUT'] = (img_name, "The image file used as input")
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

    hdu.writeto(tmp_simput_path, overwrite=True)

    # Move the last merged file to the save folder
    copy_to_save_folder(simput_file_path=tmp_simput_path, run_dir=run_dir, output_file_path=output_file_path,
                        verbose=verbose)

    # Remove the tmp simput
    # os.remove(tmp_simput_path)

    return simput_file_name
