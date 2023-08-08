import shutil
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List

from astropy.io import fits

from src.simput.utils import get_simputs
from src.sixte.image_gen import merge_ccd_eventlists, generate_fits_image, split_eventlist
from src.sixte.simulation import run_sixte_all_ccds
from utils.file_utils import decompress_gzip, compress_gzip
from utils.fits_utils import filter_evt_pattern_type
from utils.log import elog, plog
from utils.multiprocessing import mp_run
from src.xmm.pn_ff import XmmPnFf
from tempfile import TemporaryDirectory


def run_simulation(
        img_name,
        instrument_path: Path,
        simput_path: Path,
        run_dir: Path,
        res_mult,
        exposure,
        ra=0.0, dec=0.0,
        rollangle=0.0,
        debug=False, run_sim=True,
        # noise=True,
        sim_separate_ccds=False,
        verbose=True
):
    # TODO Avoid creating XmmPnFf for every single simput file
    xmm = XmmPnFf(res_mult=res_mult, sim_separate_ccds=sim_separate_ccds, verbose=verbose)
    if debug:
        xmm.pprint_xml()

    # Save the xml to the sixte instrument location
    xmm.save_xml(folder_location=instrument_path)

    # Rename the simput file to something short
    short_simput_path = run_dir / "in.simput"
    shutil.copyfile(simput_path, short_simput_path)

    if run_sim:
        run_sixte_all_ccds(xmm=xmm,
                           simput_file_path=short_simput_path,
                           rundir=run_dir,
                           exposure=exposure,
                           ra=ra,
                           dec=dec,
                           rollangle=rollangle,
                           verbose=verbose)

    # Merge all the ccd.py eventlists into one eventlist
    if sim_separate_ccds:
        combined_evt_path = merge_ccd_eventlists(input_folder=run_dir, out_filename="ccd_combined_evt.fits")
    else:
        combined_evt_path = run_dir / "ccd00_evt.fits"

    # Remove everything with pattern type (TYPE) > 4
    combined_evt_path = filter_evt_pattern_type(combined_evt_path, max_pattern_type=4, verbose=verbose)

    # split the eventlist
    split_exposure_evt_files = split_eventlist(run_dir=run_dir, eventlist_path=combined_evt_path, multiples=10000,
                                               verbose=verbose)

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
    # cdelt1 = -1.0 * xmm.cdelt_x  # / res_mult
    cdelt1 = xmm.cdelt_x  # / res_mult
    cdelt2 = xmm.cdelt_y  # / res_mult
    split_img_paths_exps = []
    for split_dict in split_exposure_evt_files:
        split_evt_file = split_dict['filename']
        split_name = split_dict['base_name']
        t_start = split_dict['t_start']
        t_stop = split_dict['t_stop']
        split_num = split_dict['split_num']
        total_splits = split_dict['total_splits']
        split_exposure = split_dict['exposure']

        final_img_name = f"{img_name}_{split_name}.fits"

        generate_fits_image(evt_file=split_evt_file, input_folder=run_dir, out_name=final_img_name, naxis1=naxis1,
                            naxis2=naxis2, crval1=crval1, crval2=crval2, crpix1=crpix1,
                            crpix2=crpix2, cdelt1=cdelt1, cdelt2=cdelt2, verbose=verbose)

        # Rename the final file, due to lengths of names we have to do this with shutil.move
        final_img_path = run_dir / final_img_name

        split_img_paths_exps.append((final_img_path, split_exposure))

        # Add specifics to the simput file
        with fits.open(final_img_path, mode="update") as hdu:
            header = hdu['PRIMARY'].header
            header['EXPOSURE'] = (split_exposure, "Exposure in seconds")
            header['ORIG_EXP'] = (exposure, "Original exposure before split in seconds")
            header['SPLIT_N'] = (split_num, "Split number (starts at 0)")
            header['SPLITS'] = (total_splits, "Total number of splits")
            header['TSTART'] = (t_start, "Start-time of split exposure")
            header['TSTOP'] = (t_stop, "Stop-time of split exposure")
            header['SIMPUT'] = (simput_path.name, "Simput used as input")
            header['RESMULT'] = (res_mult, "Resolution multiplier relative to real XMM")

            header['COMMENT'] = f"Created by Sam Sweere (samsweere@gmail.com) for ESAC at " \
                                f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}"

    return split_img_paths_exps


def run_xmm_simulation(
        simput_file: Path,
        mode,
        working_dir: Path,
        res_mult,
        exposure,
        config
):
    # TODO Path
    dataset_dir = Path("~/Documents/ESA/data/sim").expanduser() / "xmm_sim_dataset"
    debug = config["environment"]["debug"]
    verbose = config["environment"]["verbose"]

    # plog(f"Running simulation with:", verbose)
    # plog(f"Mode: {mode}", verbose)
    # plog(f"Exposure: {exposure}", verbose)
    # plog(f"Resolution multiplier: {res_mult}", verbose)
    # plog(f"Simput: {simput_path}", verbose)

    simput_filename = simput_file.name.replace(".simput.gz", "")

    res_str = f"{res_mult}x"
    img_name = f"{simput_filename}_mult_{res_mult}"

    if not debug:
        # Do not skip in debug mode
        if mode != "background":
            # Check if the max exposure time already exists, if so we know that this one was already generated
            max_exp_str = f"{round(exposure / 1000)}ks"

            # The final name always is part 0 out of 0
            final_compressed_img_name = f"{img_name}_{max_exp_str}_p_0-0.fits.gz"

            # # Save the files in a exposuretime -> mode -> resolution mult structure
            final_img_directory = dataset_dir / max_exp_str / mode / res_str

            final_compressed_file_path = final_img_directory / final_compressed_img_name
            if final_compressed_file_path.exists():
                plog(f"Image: {final_compressed_file_path} already exists, skipping", verbose=verbose)
                return final_compressed_file_path

    tmp_dir = working_dir / "xmm_sim_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory(dir=tmp_dir) as tmp:
        run_dir = Path(tmp)

        # File does not exists yet, make the run dir, unpack the simput file and run the simulation
        # Create the run_dir for this specific resolution
        # We create it here such that if the simulation fails or the file already exists
        # we have no empty run dir
        # Extract simput file
        simput_file = decompress_gzip(in_file_path=simput_file, out_file_dir=run_dir)

        # Run the simulation
        tmp_split_img_paths_exps = run_simulation(img_name=img_name,
                                                  instrument_path=Path(config["sixte"]['pn_instrument_dir']),
                                                  simput_path=simput_file,
                                                  run_dir=run_dir, res_mult=res_mult,
                                                  exposure=exposure, ra=0.0, dec=0.0, rollangle=0.0,
                                                  debug=debug, run_sim=True,
                                                  sim_separate_ccds=config["xmm"]['sim_seperate_ccds'],
                                                  verbose=verbose)

        for i, p in enumerate(tmp_split_img_paths_exps):
            file_path = p[0]
            split_exp = p[1]
            if mode == "background":
                # Remove the part numbers since they do not matter for the background
                # Split on the part numbering
                bg_filename = file_path.name
                bg_filename = f"{bg_filename.split('ks_p')[0]}ks"
                # Since we want different backgrounds we need to add an unique identifier
                # in this case the current time
                bg_filename = f"{bg_filename}_{str(time.time()).replace('.', '').ljust(17, '0')}.fits"
                new_bg_path = file_path.parent / bg_filename
                # Rename the file
                file_path.rename(new_bg_path)
                file_path = new_bg_path
            compressed_file_path = file_path.parent / f"{file_path.name}.gz"
            compress_gzip(in_file_path=file_path, out_file_path=compressed_file_path)

            if not debug:
                exp_str = str(round(split_exp / 1000)) + "ks"
                final_img_directory = dataset_dir / exp_str / mode / res_str
                final_img_directory.mkdir(parents=True, exist_ok=True)
                final_compressed_file_path = final_img_directory / compressed_file_path.name
                compressed_file_path.rename(final_compressed_file_path)
                plog(f"Saved compressed fits file at: {final_compressed_file_path}", verbose=verbose)


def run_simulation_simput(working_directory: Path, simput_file: Path, mode, config):
    for exposure in config["sixte"]['exposure']:
        # Run the simulation for all the specified resolutions
        for res_mult in config["sixte"]['res_mult']:
            try:
                # Run the simulation
                run_xmm_simulation(simput_file=simput_file, mode=mode, working_dir=working_directory,
                                   res_mult=res_mult, exposure=exposure, config=config)

            except Exception as e:
                e_m = f"ERROR: failed to run simulation with simput file: {simput_file}; exposure: {exposure}; " \
                      f"res_mult: {res_mult}; with error: {e}"
                print(e_m)

                if config["environment"]['debug']:
                    raise


def run_simulation_modes(
        mode_amount_dict: Dict[str, int],
        working_directory: Path,
        simput_path: Path,
        order: str,
        config
) -> None:
    mode_simput_files = get_simputs(simput_path=simput_path, mode_amount_dict=mode_amount_dict, order=order)

    if "background" in mode_simput_files:
        # Since we only have one background but want multiple simulations of it repeat it
        mode_simput_files["background"] = mode_simput_files["background"] * mode_amount_dict["background"]

    argument_list = []
    for mode, files in mode_simput_files.items():
        for file in files:
            argument_list.append((file, mode, config))

    # TODO Create XmmPnFf BEFORE starting the mp.Pool

    for args in argument_list:
        run_simulation_simput(working_directory, *args)

    # TODO Add multiprocessing version
