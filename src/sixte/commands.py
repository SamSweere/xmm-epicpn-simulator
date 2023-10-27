import os
from pathlib import Path
from typing import Literal, List, Union

import heasoftpy as hsp

from src.xmm_utils.external_run import run_command


def imgev(
        evt_file: Path,
        image: Path,
        coordinate_system: Literal[0, 1],
        cunit1: str,
        cunit2: str,
        naxis1: int,
        naxis2: int,
        crval1: float,
        crval2: float,
        crpix1: float,
        crpix2: float,
        cdelt1: float,
        cdelt2: float,
        clobber: bool = True,
        chatter: int = 0,
        history: bool = True
) -> None:
    params = [f"EvtFile={evt_file.resolve()}", f"Image={image.resolve()}", f"CoordinateSystem={coordinate_system}",
              f"NAXIS1={naxis1}", f"NAXIS2={naxis2}", f"CUNIT1={cunit1}", f"CUNIT2={cunit2}", f"CRVAL1={crval1}",
              f"CRVAL2={crval2}", f"CRPIX1={crpix1}", f"CRPIX2={crpix2}", f"CDELT1={cdelt1}", f"CDELT2={cdelt2}",
              f"chatter={chatter}", f"clobber={str(clobber).lower()}", f"history={str(history).lower()}"]

    cmd = f"imgev {' '.join(params)}"
    with hsp.utils.local_pfiles_context():
        run_command(cmd=cmd, cmd_input=os.devnull, verbose=True)


def runsixt(
        raw_data: Path,
        evt_file: Path,
        xml_file: Path,
        ra: float,
        dec: float,
        rollangle: float,
        simput: Path,
        exposure: int,
        chatter: int = 0,
        clobber: bool = True,
        history: bool = True
) -> None:
    params = [f"Rawdata={raw_data.resolve()}", f"EvtFile={evt_file.resolve()}", f"XMLFile={xml_file.resolve()}",
              f"RA={ra}", f"DEC={dec}", f"rollanlge={rollangle}", f"Simput={simput.resolve()}", f"Exposure={exposure}",
              f"chatter={chatter}", f"clobber={str(clobber).lower()}", f"history={str(history).lower()}"]

    cmd = f"runsixt {' '.join(params)}"
    with hsp.utils.local_pfiles_context():
        run_command(cmd=cmd, cmd_input=os.devnull, verbose=True)


def simputfile(
        simput: Path,
        ra: float = 0.0,
        dec: float = 0.0,
        src_flux: float = 0.0,
        emin: float = 1.0,
        emax: float = 10.0,
        xspec_file: Path = None,
        ascii_file: Path = None,
        image_file: Path = None,
        clobber: bool = True,
        chatter: int = 0,
        history: bool = True

) -> None:
    xspec_file = "none" if xspec_file is None else f"XSPECFile={xspec_file.resolve()}"
    ascii_file = "none" if ascii_file is None else f"ASCIIFile={ascii_file.resolve()}"
    image_file = "none" if image_file is None else f"ImageFile={image_file.resolve()}"
    params = [f"Simput={simput.resolve()}", f"RA={ra}", f"DEC={dec}", f"Emin={emin}", f"srcFlux={src_flux}",
              f"Emax={emax}", f"clobber={str(clobber).lower()}", xspec_file, ascii_file, image_file,
              f"chatter={chatter}", f"history={str(history).lower()}"]

    cmd = f"simputfile {' '.join(params)}"
    cmd = cmd.replace("none ", "")
    with hsp.utils.local_pfiles_context():
        run_command(cmd=cmd, verbose=True)


def simputmerge(
        infiles: Union[List[Path], Path],
        outfile: Path,
        fetch_extension: bool = True
) -> None:
    if isinstance(infiles, list):
        infiles = [str(f.resolve()) for f in infiles]
        infiles = ",".join(infiles)
    else:
        infiles = f"{infiles.resolve()}"
    params = [f"Infiles={infiles}", f"Outfile={outfile.resolve()}",
              f"FetchExtension={'yes' if fetch_extension else 'no'}"]

    cmd = f"simputmerge {' '.join(params)}"
    with hsp.utils.local_pfiles_context():
        run_command(cmd=cmd, cmd_input=os.devnull, verbose=True)
