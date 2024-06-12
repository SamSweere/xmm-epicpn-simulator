from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal
from uuid import uuid4

from astropy.io import fits
from loguru import logger

from src.sixte import commands
from src.sixte.image_gen import merge_ccd_eventlists, split_eventlist
from src.xmm.utils import get_cdelt, get_naxis12, get_xml_file
from src.xmm_utils.file_utils import compress_gzip, filter_event_pattern
import numpy as np 


def run_simulation(
    xml_dir: Path,
    instrument_name: Literal["epn", "emos1", "emos2"],
    xmm_filter: Literal["thin", "med", "thick"],
    simput_path: Path,
    run_dir: Path,
    res_mult: int,
    exposure: int,
    max_event_pattern: int,
    ra: float = 0.0,
    dec: float = 0.0,
    rollangle: float = 0.0,
    sim_separate_ccds: bool = False,
    consume_data: bool = True,
):
    xml_file = get_xml_file(
        xml_dir=xml_dir,
        instrument_name=instrument_name,
        res_mult=res_mult,
        xmm_filter=xmm_filter,
        sim_separate_ccds=sim_separate_ccds,
    )

    if not xml_file:
        raise FileNotFoundError(
            f"It looks like you have not created the corresponding XML file for instrument " f"'{instrument_name}'"
        )

    commands.sixtesim(
        output_path=run_dir.resolve(),
        xml_file=xml_file.resolve(),
        ra=ra,
        dec=dec,
        rollangle=rollangle,
        simput=simput_path,
        exposure=exposure,
    )

    for raw_filepath in run_dir.glob("*_raw.fits"):
        raw_filepath.unlink()

    evt_filepaths = []
    for evt_filepath in run_dir.glob("*_none"):
        evt_filepaths.append(filter_event_pattern(eventlist_path=evt_filepath, max_event_pattern=max_event_pattern))

    merged = merge_ccd_eventlists(evt_filepaths, run_dir, consume_data)

    # split the eventlist
    split_exposure_evt_files = split_eventlist(
        run_dir=run_dir,
        eventlist_path=merged,
        consume_data=consume_data,
        multiples=10000,
    )

    # See https://www.sternwarte.uni-erlangen.de/research/sixte/data/simulator_manual_v1.3.11.pdf for information
    naxis1, naxis2 = get_naxis12(instrument_name=instrument_name, res_mult=res_mult)
    cdelt1 = get_cdelt(instrument_name=instrument_name, res_mult=res_mult)
    cdelt2 = -cdelt1

    if instrument_name == "epn":
        from src.xmm.epn import get_shift_xy

        shift_y, shift_x = get_shift_xy(res_mult=res_mult)
        crpix1 = round(((naxis1 + 1) / 2.0) - shift_x, 6)
        crpix2 = round(((naxis2 + 1) / 2.0) + shift_y, 6)
    else:
        crpix1 = round(((naxis1 + 1) / 2.0), 6)
        crpix2 = round(((naxis2 + 1) / 2.0), 6)

    img_name = f"{simput_path.name.replace('.simput.gz', '')}_mult_{res_mult}"
    split_img_paths_exps = []
    for split_dict in split_exposure_evt_files:
        split_evt_file: Path = split_dict["outfile"]
        split_name = split_dict["base_name"]
        t_start = split_dict["t_start"]
        t_stop = split_dict["t_stop"]
        split_num = split_dict["split_num"]
        total_splits = split_dict["total_splits"]
        split_exposure = split_dict["exposure"]

        final_img_name = f"{img_name}_{split_name}.fits"
        final_img_path = run_dir / final_img_name

        commands.imgev(
            evt_file=split_evt_file,
            image=final_img_path,
            coordinate_system=0,
            cunit1="deg",
            cunit2="deg",
            naxis1=naxis1,
            naxis2=naxis2,
            crval1=dec,
            crval2=ra,
            crpix1=crpix1,
            crpix2=crpix2,
            cdelt1=cdelt1,
            cdelt2=cdelt2,
        )

        if consume_data:
            split_evt_file.unlink()

        split_img_paths_exps.append((final_img_path, split_exposure))

        # Add specifics to the simput file
        with fits.open(final_img_path, mode="update") as hdu:
           
            header = hdu["PRIMARY"].header
            header["EXPOSURE"] = (split_exposure, "Exposure in seconds")
            header["ORIG_EXP"] = (exposure, "Original exposure before split in seconds")
            header["SPLIT_N"] = (split_num, "Split number (starts at 0)")
            header["SPLITS"] = (total_splits, "Total number of splits")
            header["TSTART"] = (t_start, "Start-time of split exposure")
            header["TSTOP"] = (t_stop, "Stop-time of split exposure")
            header["SIMPUT"] = (simput_path.name, "Simput used as input")
            header["RESMULT"] = (res_mult, "Resolution multiplier relative to real XMM")

            header["COMMENT"] = (
                f"Created by Sam Sweere (samsweere@gmail.com) for ESAC at "
                f"{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}"
            )
            
            # Add the detector mask #JUST A TEST, REMOVE AGAIN11####
            data = hdu["PRIMARY"].data
            det_mask_path = f'res/detector_masks/pn_mask_500_2000_detxy_{res_mult}x.ds'
            with fits.open(det_mask_path) as det_mask:
                det_mask = det_mask[0].data
                
                # Compute difference in size 
                y_diff =  det_mask.shape[0] - data.shape[0]
                y_top_pad = int(np.floor(y_diff / 2.0))
                y_bottom_pad = y_diff - y_top_pad

                x_diff = det_mask.shape[1] - data.shape[1] 
                x_left_pad = int(np.floor(x_diff / 2.0))
                x_right_pad = x_diff - x_left_pad
                
             
                
                
                if y_diff >= 0:
                    # Pad the image in the y direction
                    data = np.pad(
                        data, ((y_top_pad, y_bottom_pad), (0, 0)), "constant", constant_values=0.0
                    )
                else:
                    # Crop the image in the y direction
                    data = data[abs(y_top_pad) : data.shape[0] - abs(y_bottom_pad)]

                if x_diff >= 0:
                    # Pad the image in the x direction
                    data = np.pad(
                        data, ((0, 0), (x_left_pad, x_right_pad)), "constant", constant_values=0.0
                    )
                else:
                    # Crop the image in the x direction
                    data = data[:, abs(x_left_pad) : data.shape[1] - abs(x_right_pad)]
       
                        
                
                data = data * det_mask
            
            
            hdu["PRIMARY"].data = data

    return split_img_paths_exps


def run_xmm_simulation(
    instrument_name: Literal["epn", "emos1", "emos2"],
    xml_dir: Path,
    simput_file: Path,
    mode: str,
    tmp_dir: Path,
    out_dir: Path,
    res_mult: int,
    max_event_pattern: int,
    exposure: int,
    xmm_filter: Literal["thin", "med", "thick"],
    sim_separate_ccds: bool,
    consume_data: bool,
):
    logger.info(f"Running simulations for {simput_file.resolve()}")
    with TemporaryDirectory(dir=tmp_dir) as tmp:
        run_dir = Path(tmp)
        # File does not exist yet, make the run dir, unpack the simput file and run the simulation
        # Create the run_dir for this specific resolution
        # We create it here such that if the simulation fails or the file already exists
        # we have no empty run dir

        # Run the simulation
        tmp_split_img_paths_exps = run_simulation(
            xml_dir=xml_dir,
            instrument_name=instrument_name,
            xmm_filter=xmm_filter,
            simput_path=simput_file,
            run_dir=run_dir,
            res_mult=res_mult,
            max_event_pattern=max_event_pattern,
            exposure=exposure,
            sim_separate_ccds=sim_separate_ccds,
            consume_data=consume_data,
        )

        for p in tmp_split_img_paths_exps:
            file_path: Path = p[0]
            split_exp = p[1]

            final_img_directory = out_dir / f"{round(split_exp / 1000)}ks" / mode

            if mode == "img":
                tng_name = simput_file.parts[-3]
                snapshot_num = simput_file.parts[-2]
                final_img_directory = final_img_directory / tng_name / snapshot_num

            final_img_directory = final_img_directory / f"{res_mult}x"
            final_img_directory.mkdir(parents=True, exist_ok=True)

            if mode == "bkg":
                # Remove the part numbers since they do not matter for the background
                # Split on the part numbering
                bg_filename = file_path.name
                bg_filename = f"{bg_filename.split('ks_p')[0]}ks"
                # Since we want different backgrounds we need to add an unique identifier
                bg_filename = f"{bg_filename}_{uuid4().int}.fits"
                new_bg_path = file_path.parent / bg_filename
                # Rename the file
                file_path.rename(new_bg_path)
                file_path = new_bg_path
            final_compressed_file_path = final_img_directory / f"{file_path.name}.gz"
            compress_gzip(in_file_path=file_path, out_file_path=final_compressed_file_path)
            file_path.unlink()
