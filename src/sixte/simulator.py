import time
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal, List

from astropy.io import fits
from loguru import logger

from src.sixte import commands
from src.sixte.image_gen import merge_ccd_eventlists, split_eventlist
from src.xmm.utils import get_xml_files, get_width_height_for_instrument, get_cdelt_for_instrument, \
    get_crpix12_for_instrument
from src.xmm_utils.file_utils import decompress_gzip, compress_gzip
from src.xmm_utils.fits_utils import filter_evt_pattern_type


def handle_error(error):
    logger.exception(error)


def run_simulation(
        img_name,
        instrument_name: Literal["epn", "emos1", "emos2"],
        xmm_filter: Literal["thin", "med", "thick"],
        simput_path: Path,
        run_dir: Path,
        res_mult: int,
        exposure: int,
        ra: float = 0.0,
        dec: float = 0.0,
        rollangle: float = 0.0,
        debug=False,
        sim_separate_ccds: bool = False,
        verbose: bool = True
):
    xml_paths = get_xml_files(instrument_name=instrument_name, res_mult=res_mult, xmm_filter=xmm_filter,
                              sim_separate_ccds=sim_separate_ccds)

    if not xml_paths:
        raise FileNotFoundError(f"It looks like you have not created the corresponding XML files for instrument "
                                f"'{instrument_name}'")

    evt_filepaths: List[Path] = []
    for xml_path in xml_paths:
        ccd_name = xml_path.stem
        evt_filepath = run_dir / f"{ccd_name}_evt.fits"
        commands.runsixt(raw_data=run_dir / f"{ccd_name}_raw.fits", evt_file=evt_filepath, xml_file=xml_path.resolve(),
                         ra=ra, dec=dec, rollangle=rollangle, simput=simput_path, exposure=exposure)
        evt_filepaths.append(evt_filepath)

    if len(evt_filepaths) > 1:
        # Merge all the ccd.py eventlists into one eventlist
        evt_filepath = merge_ccd_eventlists(infiles=evt_filepaths, out_dir=run_dir, verbose=verbose)
    else:
        evt_filepath = evt_filepaths[0]

    # Remove everything with pattern type (TYPE) > 4
    combined_evt_path = filter_evt_pattern_type(evt_filepath, max_pattern_type=4, verbose=verbose)

    # split the eventlist
    split_exposure_evt_files = split_eventlist(run_dir=run_dir, eventlist_path=combined_evt_path, multiples=10000,
                                               verbose=verbose)

    # See https://www.sternwarte.uni-erlangen.de/research/sixte/data/simulator_manual_v1.3.11.pdf for information
    naxis2, naxis1 = get_width_height_for_instrument(instrument_name=instrument_name, res_mult=res_mult)
    cdelt1 = cdelt2 = get_cdelt_for_instrument(instrument_name=instrument_name, res_mult=res_mult)
    crpix2, crpix1 = get_crpix12_for_instrument(instrument_name=instrument_name, res_mult=res_mult)

    split_img_paths_exps = []
    for split_dict in split_exposure_evt_files:
        split_evt_file = split_dict['outfile']
        split_name = split_dict['base_name']
        t_start = split_dict['t_start']
        t_stop = split_dict['t_stop']
        split_num = split_dict['split_num']
        total_splits = split_dict['total_splits']
        split_exposure = split_dict['exposure']

        final_img_name = f"{img_name}_{split_name}.fits"
        final_img_path = run_dir / final_img_name

        commands.imgev(evt_file=split_evt_file, image=final_img_path, coordinate_system=0, cunit1="deg", cunit2="deg",
                       naxis1=naxis1, naxis2=naxis2, crval1=dec, crval2=ra, crpix1=crpix1, crpix2=crpix2, cdelt1=cdelt1,
                       cdelt2=cdelt2)

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
        img_name: str,
        mode: str,
        tmp_dir: Path,
        out_dir: Path,
        res_mult: int,
        exposure: int,
        instrument_name: Literal["epn", "emos1", "emos2"],
        xmm_filter: Literal["thin", "med", "thick"],
        sim_separate_ccds: bool,
        debug: bool = False,
        verbose: bool = True
):
    logger.info(f"Running simulations for {simput_file.resolve()}")
    with TemporaryDirectory(dir=tmp_dir) as tmp:
        run_dir = Path(tmp)
        # File does not exist yet, make the run dir, unpack the simput file and run the simulation
        # Create the run_dir for this specific resolution
        # We create it here such that if the simulation fails or the file already exists
        # we have no empty run dir
        # Extract simput file
        simput_file = decompress_gzip(in_file_path=simput_file, out_file_dir=run_dir)

        # Run the simulation
        tmp_split_img_paths_exps = run_simulation(img_name=img_name,
                                                  instrument_name=instrument_name,
                                                  xmm_filter=xmm_filter,
                                                  simput_path=simput_file,
                                                  run_dir=run_dir, res_mult=res_mult,
                                                  exposure=exposure,
                                                  debug=debug,
                                                  sim_separate_ccds=sim_separate_ccds,
                                                  verbose=verbose)

        for i, p in enumerate(tmp_split_img_paths_exps):
            file_path = p[0]
            split_exp = p[1]

            final_img_directory = out_dir / f"{round(split_exp / 1000)}ks" / mode / f"{res_mult}x"
            final_img_directory.mkdir(parents=True, exist_ok=True)

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
            final_compressed_file_path = final_img_directory / f"{file_path.name}.gz"
            compress_gzip(in_file_path=file_path, out_file_path=final_compressed_file_path)
