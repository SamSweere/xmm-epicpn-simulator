import os
import shutil
import time
from datetime import datetime

from astropy.io import fits
from tqdm import tqdm

from simput.utils import get_simputs
from sixte.image_gen import merge_ccd_eventlists, generate_fits_image, split_eventlist
from sixte.simulation import run_sixte_all_ccds
from utils.file_utils import decompress_gzip, compress_gzip
from utils.fits_utils import filter_evt_pattern_type
from utils.log import elog, plog
from utils.multiprocessing import get_multiprocessing_pool
from xmm.pn_ff import XmmPnFf


def run_simulation(img_name, instrument_path, simput_path, run_dir, res_mult, exposure,
                   ra=0.0, dec=0.0,
                   rollangle=0.0,
                   debug=False, run_sim=True, noise=True, sim_separate_ccds=False, verbose=True):
    xmm = XmmPnFf(res_mult=res_mult, noise=noise, sim_separate_ccds=sim_separate_ccds, verbose=verbose)
    if debug:
        xmm.pprint_xml()

    # Save the xml to the sixte instrument location
    xmm.save_xml(folder_location=instrument_path)

    if not os.path.isdir(run_dir):
        e_message = f'The output root folder {run_dir} does not exist. Cannot continue.'
        elog(e_message)
        raise FileNotFoundError(e_message)

    # Rename the simput file to something short
    short_simput_name = "in.simput"
    short_simput_path = os.path.join(run_dir, short_simput_name)
    shutil.copyfile(simput_path, short_simput_path)

    if run_sim:
        run_sixte_all_ccds(xmm=xmm, simput_file_path=short_simput_path, rundir=run_dir, exposure=exposure, ra=ra,
                           dec=dec,
                           rollangle=rollangle, verbose=verbose)

    # Merge all the ccd eventlists into one eventlist
    if sim_separate_ccds:
        combined_evt_path = merge_ccd_eventlists(input_folder=run_dir)
    else:
        combined_evt_path = os.path.join(run_dir, "ccd00_evt.fits")

    # Remove everything with pattern type (TYPE) > 4
    combined_evt_path = filter_evt_pattern_type(combined_evt_path, max_pattern_type=4, verbose=verbose)

    # split the eventlist
    split_exposure_evt_files = split_eventlist(run_dir=run_dir, eventlist_path=combined_evt_path, multiples=10000,
                                               verbose=verbose)

    split_img_paths_exps = []

    for split_dict in split_exposure_evt_files:
        split_evt_file = split_dict['filename']
        split_name = split_dict['base_name']
        t_start = split_dict['t_start']
        t_stop = split_dict['t_stop']
        split_num = split_dict['split_num']
        total_splits = split_dict['total_splits']
        split_exposure = split_dict['exposure']

        final_img_name = img_name + "_" + split_name + ".fits"

        # See https://www.sternwarte.uni-erlangen.de/research/sixte/data/simulator_manual_v1.3.11.pdf for information
        # NAXIS1 and NAXIS2 keywords give the number of x- andy-pixels in the calculated exposure map
        # CRVAL1 and CRVAL2 define the center of the map
        # CRPIX1 and CRPIX2 are the point corresponding to the optical axis on the instrument define the corresponding
        # pixel coordinate
        # CDELT1 and CDELT2 give the pixel sizes in degrees
        naxis2 = round(xmm.total_pixel_height)  # * res_mult
        naxis1 = round(xmm.total_pixel_width)  # * res_mult
        crval1 = dec
        crval2 = ra
        crpix1 = round(xmm.crpix1, 12)  # * res_mult
        crpix2 = round(xmm.crpix2, 12)  # * res_mult
        cdelt1 = -1.0 * xmm.cdelt_x  # / res_mult
        cdelt2 = xmm.cdelt_y  # / res_mult
        # final_img_name = f"{run_id}.fits"
        run_img_name = "img_comb.fits"

        generate_fits_image(evt_file=split_evt_file, input_folder=run_dir, out_name=run_img_name, naxis1=naxis1,
                            naxis2=naxis2, crval1=crval1, crval2=crval2, crpix1=crpix1,
                            crpix2=crpix2, cdelt1=cdelt1, cdelt2=cdelt2, verbose=verbose)

        # Rename the final file, due to lengths of names we have to do this with shutil.move
        final_img_path = os.path.join(run_dir, final_img_name)
        shutil.move(os.path.join(run_dir, run_img_name), final_img_path)

        # Add specifics to the simput file
        hdu = fits.open(final_img_path)
        header = hdu['PRIMARY'].header
        header['EXPOSURE'] = (split_exposure, "Exposure in seconds")
        header['ORIG_EXP'] = (exposure, "Original exposure before split in seconds")
        header['SPLIT_N'] = (split_num, "Split number (starts at 0)")
        header['SPLITS'] = (total_splits, "Total number of splits")
        header['TSTART'] = (t_start, "Start-time of split exposure")
        header['TSTOP'] = (t_stop, "Stop-time of split exposure")
        header['SIMPUT'] = (os.path.basename(simput_path), "Simput used as input")
        header['RESMULT'] = (res_mult, "Resolution multiplier relative to real XMM")

        header['COMMENT'] = f"Created by Sam Sweere (samsweere@gmail.com) for ESAC at " \
                            f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}"

        hdu.writeto(final_img_path, overwrite=True)

        split_img_paths_exps.append((final_img_path, split_exposure))

    return split_img_paths_exps


def run_xmm_simulation(simput_path, mode, run_dir, res_mult, exposure, config):
    plog(f"Running simulation with:", config['verbose'])
    plog(f"Mode: {mode}", config['verbose'])
    plog(f"Exposure: {exposure}", config['verbose'])
    plog(f"Resolution multiplier: {res_mult}", config['verbose'])
    plog(f"Simput: {simput_path}", config['verbose'])

    simput_filename = os.path.basename(simput_path).replace(".simput.gz", "")

    res_str = str(res_mult) + "x"

    img_name = simput_filename + "_mult_" + str(res_mult)

    if not config['debug']:
        # Do not skip in debug mode
        if mode != "background":
            # Check if the max exposure time already exists, if so we know that this one was already generated
            max_exp_str = str(round(exposure / 1000)) + "ks"

            # The final name always is part 0 out of 0
            final_compressed_img_name = img_name + f"_{max_exp_str}_p_0-0" + ".fits.gz"
            # compressed_file_path = os.path.join(run_dir, final_compressed_img_name)

            # # Save the files in a exposuretime -> mode -> resolution mult structure
            final_img_directory = os.path.join(config['dataset_dir'], max_exp_str, mode, res_str)

            final_compressed_file_path = os.path.join(final_img_directory, final_compressed_img_name)
            if os.path.exists(final_compressed_file_path):
                plog(f"Image: {final_compressed_file_path} already exists, skipping", verbose=config['verbose'])
                return final_compressed_file_path

    # File does not exists yet, make the run dir, unpack the simput file and run the simulation
    # Create the run_dir for this specific resolution
    # We create it here such that if the simulation fails or the file already exists
    # we have no empty run dir
    time_str = str(time.time()).replace(".", "")
    run_dir = os.path.join(run_dir, str(time_str))
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
    else:
        m = f"Error, rundir {run_dir} already exists"
        elog(m)
        raise FileExistsError(m)

    # Extract simput file
    simput_path = decompress_gzip(in_file_path=simput_path, out_file_dir=run_dir)

    # Run the simulation
    tmp_split_img_paths_exps = run_simulation(img_name=img_name,
                                              instrument_path=config['instrument_dir'], simput_path=simput_path,
                                              run_dir=run_dir, res_mult=res_mult,
                                              exposure=exposure, ra=0.0, dec=0.0, rollangle=0.0,
                                              debug=config['debug'], run_sim=True, noise=False,
                                              sim_separate_ccds=config['simulate_separate_ccds'],
                                              verbose=config['verbose'])

    if mode == "background":
        # Remove the part numbers since they do not matter for the background
        for i in range(len(tmp_split_img_paths_exps)):
            back_path = tmp_split_img_paths_exps[i][0]
            split_exp = tmp_split_img_paths_exps[i][1]

            # SPlit on the part numbering
            new_back_path = "".join(back_path.split("ks_p")[:-1]) + "ks"
            # new_back_path = "background_" + "mult" + "_" + str(res_mult) + "_" + str(round(split_exp / 1000)) + "ks"
            # Since we want different backgrounds we need to add an unique identifier
            # in this case the current time
            new_back_path += "_" + str(time.time()).replace(".", "").ljust(17, '0')
            new_back_path += ".fits"

            new_back_path = os.path.join(run_dir, new_back_path)
            # Rename the file
            shutil.move(back_path, new_back_path)
            tmp_split_img_paths_exps[i] = (new_back_path, split_exp)

    for tmp_split_img_path, split_exp in tmp_split_img_paths_exps:
        compressed_file_path = tmp_split_img_path.replace('.fits', '.fits.gz')
        exp_str = str(round(split_exp / 1000)) + "ks"
        final_img_directory = os.path.join(config['dataset_dir'], exp_str, mode, res_str)

        final_compressed_file_path = os.path.join(final_img_directory, os.path.basename(compressed_file_path))

        # Compress the fits file and move it too the correct output dir
        compress_gzip(in_file_path=tmp_split_img_path, out_file_path=compressed_file_path)

        if not config['debug']:
            # Create the final image directory (if in debug mode or if the process fails or is interrupted this
            # will thus not lead to empty directories)
            if not os.path.exists(final_img_directory):
                os.makedirs(final_img_directory)

            shutil.move(compressed_file_path, final_compressed_file_path)

            plog(f"Saved compressed fits file at: {final_compressed_file_path}", verbose=config['verbose'])

    if not config['debug']:
        # Remove the tmp rundir if it exists
        plog(f'Removing run directory {run_dir}', verbose=config['verbose'])
        # Remove the run directory
        shutil.rmtree(run_dir)

    return


def run_simulation_simput(simput_path, mode, run_dir, config):
    # Run the simulation for all the specified exposures
    for exposure in config['exposure']:
        # Run the simulation for all the specified resolutions
        for res_mult in config['res_mult']:
            try:
                # Run the simulation
                run_xmm_simulation(simput_path=simput_path, mode=mode, run_dir=run_dir,
                                   res_mult=res_mult, exposure=exposure, config=config)

            except Exception as e:
                e_m = f"ERROR: failed to run simulation with simput file: {simput_path}; exposure: {exposure}; res_mult: {res_mult}; with error: {e}"
                print(e_m)
                elog(e_m)

                if config['debug']:
                    raise


def run_simulation_mode(mode, amount, run_dir, config):
    if amount == 0:
        return

    # Get the pre-generated simputs
    simputs = get_simputs(simput_base_path=config['simput_base_path'], mode=mode, amount=amount,
                          order=config['simput_order'])

    if mode == "background":
        # Since we only have one background but want multiple simulations of it repeat it
        simputs = simputs * amount

    # Run the processes parallel with a progress bar
    pool = get_multiprocessing_pool(num_processes=config['num_processes'], gb_per_process=config['gb_per_process'])
    jobs = []

    for simput_path in simputs:
        # Run the simulation for the simput file
        job = pool.apply_async(run_simulation_simput,
                               args=(simput_path, mode, run_dir, config))
        jobs.append(job)

        # For debugging:
        # run_simulation_simput(simput_path, mode, run_dir, config)

    pool.close()
    result_list_tqdm = []
    for job in tqdm(jobs):
        result_list_tqdm.append(job.get())
