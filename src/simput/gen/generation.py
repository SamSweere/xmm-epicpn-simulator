import os
import random
import shutil
import time
from pathlib import Path

import numpy as np
from tqdm import tqdm

from simput.img_simputgen import create_img_simput
from src.simput.utils import get_spectrumfile
from src.simput.agn import get_fluxes
from src.simput.gen import generate_background, generate_exposure_map, generate_point_source
from src.simput.utils import merge_and_save_simputs
from utils.file_utils import compress_gzip
from utils.log import elog, plog


def create_background(run_dir, verbose):
    # TODO Move to /res
    spectrum_ds_file = os.path.join(os.path.dirname(__file__), "../../../simput/spectrum/pntffg_spectrum.ds")
    # spectrum_ds_file = os.path.join(os.path.dirname(__file__), "spectrum/back_spec.fits")

    emin = 0.15
    emax = 15.0

    simput_file_name = generate_background(run_dir=run_dir, spectrum_file=spectrum_ds_file, emin=emin, emax=emax,
                                           verbose=verbose)

    return simput_file_name


def create_exposure_map(run_dir, verbose):
    # Take the background spectrum as template
    # TODO Move to /res
    spectrum_ds_file = os.path.join(os.path.dirname(__file__), "../../../simput/spectrum/pntffg_spectrum.ds")

    emin = 0.15
    emax = 15.0

    simput_file_name = generate_exposure_map(run_dir=run_dir, spectrum_file=spectrum_ds_file, emin=emin,
                                             emax=emax,
                                             verbose=verbose)

    return simput_file_name


def create_random_sources(
        run_dir: Path,
        num_sources=10
):
    simput_files = []
    center = (0, 0)

    spectrum_file = get_spectrumfile(run_dir=run_dir, norm=0.01)

    # emin = 0.15, emax = 15.0 is ideal
    # For some reason if I go into the emin=0.5 and emax=2.0 sixte will not complete anymore
    emin = 0.15
    emax = 15.0

    name = f"random_n_{num_sources}_{emin}eV_p1_{emax}eV"
    output_file = run_dir / f"{name}.simput"

    # Create the point sources
    for i in range(num_sources):
        print(f"Generating point source {i}/{num_sources - 1}")
        simput_file = run_dir / f"ps_{i}.simput"
        generate_point_source(emin=emin, emax=emax, center_point=center, xspec_file=spectrum_file,
                              simput_file_path=simput_file, src_flux='random', offset='random')
        simput_files.append(simput_file)

    return merge_and_save_simputs(run_dir=run_dir, simput_files=simput_files, output_file=output_file)


def create_agn_sources(
        run_dir: Path,
        verbose=True
):
    simput_files = []
    emin = 0.5
    emax = 2.0

    # Use the current time as id, such that clashes don't happen
    unique_id = str(time.time()).replace(".", "").ljust(17, '0')
    name = f"agn_{unique_id}_p0_{emin}ev_p1_{emax}ev"
    simput_file_name = f"{name}.simput"
    output_file_path = run_dir / simput_file_name

    # Get the fluxes from the agn distribution
    fluxes = get_fluxes()

    # Get the spectrum file
    spectrum_file = get_spectrumfile(run_dir=run_dir, norm=0.001, verbose=verbose)

    for i in tqdm(range(len(fluxes)), disable=not verbose):
        flux = fluxes[i]
        if verbose:
            print("Creating source with flux:", flux)
        simput_file_path = run_dir / f"ps_{i}.simput"
        generate_point_source(emin=emin, emax=emax, simput_file_path=simput_file_path, src_flux=flux,
                              xspec_file=spectrum_file, offset='random', verbose=verbose)
        simput_files.append(simput_file_path)

    return merge_and_save_simputs(run_dir=run_dir, simput_files=simput_files, output_file=output_file_path,
                                  verbose=verbose)


# def create_img_simput(run_dir, img_path, output_file_dir, verbose=True):
#     img_name = os.path.split(img_path)[-1]
#
#     simput_files = []
#     emin = 0.5
#     emax = 2.0
#     # TODO: random fluctuate norm
#     norm = 0.01
#
#     # Use the current time as id, such that clashes don't happen
#     id = str(time.time()).replace(".", "")
#     fits_name_path = img_name.replace(".fits", "")
#     name = f"{fits_name_path}_p0_{emin}ev_p1_{emax}ev_norm_{norm}_{id}"
#     simput_file_name = name + ".simput"
#     output_file_path = os.path.join(output_file_dir, simput_file_name)
#
#     # Get the spectrum file
#     spectrum_file = get_spectrumfile(run_dir=run_dir, norm=norm, verbose=verbose)
#
#     # Epic pn has a fov of 30 arcmin = 0.5 degrees.
#     pn_fov = 0.5
#
#     # Randomly position the point source within the fov
#     center_point = (0.0, 0.0)
#
#     offset = (0.0, 0.0)
#     #TODO: randomize offset
#     # if offset == 'random':
#     #     offset = np.random.uniform(low=-1.0 * pn_fov / 2, high=pn_fov / 2, size=2)
#
#     location = (center_point[0] + offset[0], center_point[1] + offset[1])
#     ra = location[0]
#     if ra < 0:
#         ra = 360 + ra
#     dec = location[1]
#     # if dec < 0:
#     #     dec = 90 + dec
#
#     # Create a tmp simput name such that the simputfile does not crash
#     tmp_simput_path = os.path.join(run_dir, "img.simput")
#
#     # We need the xspec from headas
#     simput_command = f"simputfile RA={ra} Dec={dec} XSPECFile={spectrum_file} Emin={emin} Emax={emax} " \
#                      f"Simput={tmp_simput_path} ImageFile={img_path}"
#     run_headas_command(simput_command, verbose=verbose)
#
#     # TODO: rotate simput by a random factor
#
#     # Move the last merged file to the save folder
#     copy_to_save_folder(simput_file_path=tmp_simput_path, run_dir=run_dir, output_file_path=output_file_path,
#                         verbose=verbose)
#
#     return output_file_path


def create_test_grid(run_dir, flux=1.0e-13, step_size=10, center_offset=(0.005, 0.015), verbose=True):
    name = f"test_grid_f_{flux}_step_{step_size}"
    simput_file_name = f"{name}.simput"
    output_file_path = run_dir / simput_file_name

    simput_files = []

    spectrum_file = get_spectrumfile(run_dir=run_dir, norm=0.01, verbose=verbose)

    # emin = 0.15, emax = 15.0 is ideal
    # For some reason if I go into the emin=0.5 and emax=2.0 sixte will not complete anymore
    emin = 0.15
    emax = 15.0

    xmm_fov = 0.5
    x_loc = np.linspace(-xmm_fov / 2, xmm_fov / 2, step_size)
    y_loc = np.linspace(-xmm_fov / 2, xmm_fov / 2, step_size)
    counter = 0
    # Create the point sources
    for x in x_loc:
        for y in y_loc:
            if verbose:
                print(f"Generating point source {counter}/{step_size ** 2}")
            ps_simput_file_name = f"ps_{counter}.simput"
            simput_file_path = os.path.join(run_dir, ps_simput_file_name)
            generate_point_source(emin=emin, emax=emax, center_point=center_offset, xspec_file=spectrum_file,
                                  simput_file_path=simput_file_path, src_flux=flux, offset=(x, y), verbose=verbose)
            simput_files.append(simput_file_path)
            counter += 1

    # Create one point source at (0, 0)
    if verbose:
        print(f"Generating point source at (0, 0)")
    simput_file_path = run_dir / f"ps_center.simput"
    generate_point_source(emin=emin, emax=emax, center_point=(0.0, 0.0), xspec_file=spectrum_file,
                          simput_file_path=simput_file_path, src_flux=flux, offset=(0.0, 0.0), verbose=verbose)
    simput_files.append(simput_file_path)

    # print(f"Generating point source at (0.125, 0.125) deg")
    # simput_file_name = f"ps_extra.simput"
    # simput_file_path = os.path.join(run_dir, simput_file_name)
    # create_point_source(emin=emin, emax=emax, center_point=(0.0, 0.0), xspec_file=spectrum_file,
    #                     simput_file_path=simput_file_path, src_flux=flux, offset=(0.125, 0.125))
    # simput_files.append(simput_file_path)
    #
    merge_and_save_simputs(run_dir=run_dir, simput_files=simput_files, output_file=output_file_path,
                           verbose=verbose)

    return simput_file_name


def generate_simput(mode, img_settings, run_dir, output_file_dir, simput_config, debug, verbose):
    simput_dir = simput_config['simput_dir']
    # Create the temporary run_dir
    # Making the run_id too long will result in failed Headas commands
    run_id = str(random.randint(100000000000, 999999999999))  # str(time.time()).replace(".", "").ljust(17, '0')[-12:] +
    run_dir = os.path.join(run_dir, str(run_id))

    # output_file_name = str(run_id) + ".simput"
    # output_file_path = os.path.join(simput_dir, output_file_name)
    # # Check if file already exists
    # if check_file_exist(output_file_path):
    #     # File already exists
    #     return output_file_path

    # Simput file does not exist yet, create it
    os.makedirs(run_dir)

    if mode == "test_grid":
        # Returns simput file path
        simput_file_name = create_test_grid(run_dir=run_dir, verbose=verbose)
    elif mode == "agn":
        simput_file_name = create_agn_sources(run_dir=run_dir, verbose=verbose)
    elif mode == "img":
        simput_file_name = create_img_simput(run_dir=run_dir, img_settings=img_settings, verbose=verbose)
    elif mode == "background":
        simput_file_name = create_background(run_dir=run_dir, verbose=verbose)
    elif mode == "exposure_map":
        simput_file_name = create_exposure_map(run_dir=run_dir, verbose=verbose)
    elif mode == "random":
        simput_file_name = create_random_sources(run_dir=run_dir)
    else:
        e_message = f"Mode {mode} is not supported"
        elog(e_message)
        raise ValueError(e_message)

    # Compress the simput file and move it too the correct output dir
    tmp_simput_path = os.path.join(run_dir, simput_file_name)
    tmp_compressed_path = tmp_simput_path + ".gz"
    final_simput_path = os.path.join(output_file_dir, simput_file_name)
    compressed_file_path = final_simput_path + ".gz"
    if os.path.exists(compressed_file_path):
        m = f"WARNING: Simput file {compressed_file_path} already exists, skipping"
        plog(m, verbose)
        return compressed_file_path
    else:
        compress_gzip(in_file_path=tmp_simput_path, out_file_path=tmp_compressed_path)
        # Once compressed move the file, this will ensure that the dataset will always have the full file
        shutil.move(tmp_compressed_path, compressed_file_path)

    plog(f"Saved compressed simput file at: {compressed_file_path}", verbose)

    if not debug:
        # Remove the temporary run_dir
        shutil.rmtree(run_dir)
        plog(f"Removed tmp dir {run_dir}", verbose=verbose)

    return final_simput_path
