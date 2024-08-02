from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal
from uuid import uuid4

from loguru import logger

from src.heasoft import heasoft as hsp
from src.sixte import commands
from src.sixte.image_gen import split_eventlist
from src.xmm.tools import get_cdelt, get_crpix12, get_naxis12, get_xml_files


def run_simulation(
    tmp_dir: Path,
    out_dir: Path,
    xml_dir: Path,
    instrument_name: Literal["epn", "emos1", "emos2"],
    xmm_filter: Literal["thin", "med", "thick"],
    simput_path: Path,
    res_mults: list[int],
    exposure: int,
    max_event_pattern: int,
    mode: str,
    ra: float = 0.0,
    dec: float = 0.0,
    rollangle: float = 0.0,
    sim_separate_ccds: bool = False,
    consume_data: bool = True,
    emasks: dict[int, Path] = None,
) -> list[Path] | None:
    split_img_paths_exps = []
    logger.info(f"Running simulations for {simput_path}")
    for res_mult in res_mults:
        with TemporaryDirectory(dir=tmp_dir) as tmp:
            run_dir = Path(tmp)
            xml_files = get_xml_files(
                xml_dir=xml_dir,
                instrument_name=instrument_name,
                res_mult=res_mult,
                xmm_filter=xmm_filter,
                sim_separate_ccds=sim_separate_ccds,
            )

            logger.debug(f"Got XML files: {xml_files}")

            commands.runsixt_sixtesim(
                output_path=run_dir,
                xml_files=xml_files,
                ra=ra,
                dec=dec,
                rollangle=rollangle,
                simput=simput_path,
                exposure=exposure,
            )

            evt_filepaths = list(run_dir.glob("chip*_none"))

            assert len(evt_filepaths) > 0

            merged = hsp.ftmerge(evt_filepaths, run_dir / "merged_events.fits", consume_data)

            if max_event_pattern < 12:
                merged = hsp.ftcopy(f"{merged}[EVENTS][TYPE <= {max_event_pattern}]", merged)
                for i in range(max_event_pattern + 1, 13):
                    hsp.fthedit(f"{merged}[EVENTS]", f"NGRAD{i}", "add", "0")
                    hsp.fthedit(f"{merged}[EVENTS]", f"NPGRA{i}", "add", "0")

            # split the eventlist
            split_events = split_eventlist(
                run_dir=run_dir,
                eventlist_path=merged,
                consume_data=consume_data,
                multiples=10000,
            )

            # See https://www.sternwarte.uni-erlangen.de/research/sixte/data/simulator_manual_v1.3.11.pdf
            # for information
            naxis1, naxis2 = get_naxis12(instrument_name=instrument_name, res_mult=res_mult)
            cdelt1, cdelt2 = get_cdelt(instrument_name=instrument_name, res_mult=res_mult)
            crpix1, crpix2 = get_crpix12(instrument_name, res_mult)

            logger.debug(
                f"NAXIS1: {naxis1}\tNAXIS2: {naxis2}\tCDELT1: {cdelt1}\tCDELT2: {cdelt2}\tCRPIX1: {crpix1}"
                f"\tCRPIX2: {crpix2}"
            )

            img_name = f"{simput_path.name.replace('.simput.gz', '')}_mult_{res_mult}"

            for split_event, split_exp in split_events:
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

                logger.debug(f"Created image {final_img_path} with exposure {split_exp}")

                if consume_data:
                    split_event.unlink()

                final_img_directory = out_dir / f"{round(split_exp / 1000)}ks"

                if mode == "img":
                    tng_name = simput_path.parts[-3]
                    snapshot_num = simput_path.parts[-2]
                    final_img_directory = final_img_directory / tng_name / snapshot_num

                final_img_directory = final_img_directory / f"{res_mult}x"
                final_img_directory.mkdir(parents=True, exist_ok=True)

                if mode == "bkg":
                    # Remove the part numbers since they do not matter for the background
                    # Split on the part numbering
                    bg_filename = final_img_path.name
                    bg_filename = f"{bg_filename.split('ks_p')[0]}ks"
                    # Since we want different backgrounds we need to add an unique identifier
                    bg_filename = f"{bg_filename}_{uuid4().int}.fits"
                    new_bg_path = final_img_path.parent / bg_filename
                    # Rename the file
                    final_img_path.rename(new_bg_path)
                    final_img_path = new_bg_path
                final_compressed_file_path = final_img_directory / f"{final_img_path.name}.gz"

                # Add specifics to the simput file and apply emask if requested
                emask = emasks[res_mult]
                if emask is not None:
                    hsp.ftimgcalc(final_compressed_file_path, "A * B", a=final_img_path, b=emask)
                else:
                    hsp.ftcopy(final_img_path, final_compressed_file_path)
                final_img_path.unlink()
                split_img_paths_exps.append(final_compressed_file_path)
    return split_img_paths_exps
