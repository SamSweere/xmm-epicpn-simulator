from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal
from uuid import uuid4

import numpy as np
from astropy.io import fits
from loguru import logger

from src.sixte import commands
from src.sixte.image_gen import merge_ccd_eventlists, split_eventlist
from src.tools.files import compress_gzip, filter_event_pattern
from src.xmm.tools import get_cdelt, get_crpix12, get_naxis12, get_xml_file


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
    emask: Path = None,
) -> list[tuple[Path, int]] | None:
    xml_file = get_xml_file(
        xml_dir=xml_dir,
        instrument_name=instrument_name,
        res_mult=res_mult,
        xmm_filter=xmm_filter,
        sim_separate_ccds=sim_separate_ccds,
    )

    commands.sixtesim(
        output_path=run_dir,
        xml_file=xml_file,
        ra=ra,
        dec=dec,
        rollangle=rollangle,
        simput=simput_path,
        exposure=exposure,
    )

    for raw_filepath in run_dir.glob("*_raw.fits"):
        raw_filepath.unlink()

    logger.debug("All *_raw.fits have been deleted")

    evt_filepaths = []
    for evt_filepath in run_dir.glob("*_none"):
        evt_filepath = filter_event_pattern(eventlist_path=evt_filepath, max_event_pattern=max_event_pattern)
        if evt_filepath is not None:
            evt_filepaths.append(evt_filepath)

    if not evt_filepaths:
        logger.warning(f"No events have been detected for detector {instrument_name} for {simput_path}!")
        return None

    merged = merge_ccd_eventlists(evt_filepaths, run_dir, consume_data)

    # split the eventlist
    split_events = split_eventlist(
        run_dir=run_dir,
        eventlist_path=merged,
        consume_data=consume_data,
        multiples=10000,
    )

    # See https://www.sternwarte.uni-erlangen.de/research/sixte/data/simulator_manual_v1.3.11.pdf for information
    naxis1, naxis2 = get_naxis12(instrument_name=instrument_name, res_mult=res_mult)
    cdelt1, cdelt2 = get_cdelt(instrument_name=instrument_name, res_mult=res_mult)
    crpix1, crpix2 = get_crpix12(instrument_name, res_mult)

    logger.debug(
        f"NAXIS1: {naxis1}\tNAXIS2: {naxis2}\tCDELT1: {cdelt1}\tCDELT2: {cdelt2}\tCRPIX1: {crpix1}\tCRPIX2: {crpix2}"
    )

    img_name = f"{simput_path.name.replace('.simput.gz', '')}_mult_{res_mult}"
    if emask is not None:
        logger.info("A mask will be applied")

        with fits.open(emask, mode="readonly") as f:
            emask = f["mask"].data if "mask" in f else f[0].data

        if instrument_name == "emos1":
            emask = np.rot90(emask)

    split_img_paths_exps = []
    for split_event in split_events:
        exposure = fits.getheader(split_event, "EVENTS")["EXPOSURE"]
        final_img_path = run_dir / f"{img_name}_{split_event.stem}.fits"

        commands.imgev(
            evt_file=split_event,
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
            split_event.unlink()

        split_img_paths_exps.append((final_img_path, exposure))

        # Add specifics to the simput file and apply emask if requested
        if emask is not None:
            with fits.open(final_img_path, mode="update") as hdu:
                hdu["PRIMARY"].data = hdu["PRIMARY"].data * emask

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
    emask: Path = None,
) -> list[Path]:
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
            emask=emask,
        )

        res = []
        if tmp_split_img_paths_exps is None:
            logger.warning(f"Something went wrong with {simput_file}")
            return res

        for p in tmp_split_img_paths_exps:
            file_path: Path = p[0]
            split_exp = p[1]

            final_img_directory = out_dir / mode / f"{round(split_exp / 1000)}ks"

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
            compress_gzip(in_file_path=file_path, out_file_path=final_compressed_file_path, remove_file=True)
            res.append(final_compressed_file_path)
    return res
