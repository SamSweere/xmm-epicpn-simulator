from pathlib import Path
from typing import Literal

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
    history: bool = True,
) -> None:
    imgev = (
        f"imgev EvtFile={evt_file.resolve()} Image={image.resolve()} CoordinateSystem={coordinate_system} "
        f"NAXIS1={naxis1} NAXIS2={naxis2} CUNIT1={cunit1} CUNIT2={cunit2} CRVAL1={crval1} CRVAL2={crval2} "
        f"CRPIX1={crpix1} CRPIX2={crpix2} CDELT1={cdelt1} CDELT2={cdelt2} chatter={chatter} "
        f"clobber={str(clobber).lower()} history={str(history).lower()}"
    )
    run_command(cmd=imgev)


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
    history: bool = True,
) -> None:
    runsixt = (
        f"runsixt Rawdata={raw_data.resolve()} EvtFile={evt_file.resolve()} XMLFile={xml_file.resolve()} "
        f"RA={ra} DEC={dec} rollangle={rollangle} Simput={simput.resolve()} Exposure={exposure} "
        f"chatter={chatter} clobber={str(clobber).lower()} history={str(history).lower()}"
    )
    run_command(cmd=runsixt)


def simputfile(
    simput: Path,
    ra: float = 0.0,
    dec: float = 0.0,
    src_flux: float = 0.0,
    emin: float = 1.0,
    emax: float = 10.0,
    xspec_file: Path | None = None,
    ascii_file: Path | None = None,
    image_file: Path | None = None,
    clobber: bool = True,
    chatter: int = 0,
    history: bool = True,
) -> None:
    simputfile = (
        f"simputfile Simput={simput.resolve()} RA={ra} DEC={dec} srcFlux={src_flux} Emin={emin} "
        f"Emax={emax} clobber={str(clobber).lower()} chatter={chatter} "
        f"history={str(history).lower()}"
    )

    if xspec_file is not None:
        simputfile += f" XSPECFile={xspec_file.resolve()}"

    if ascii_file is not None:
        simputfile += f" ASCIIFile={ascii_file.resolve()}"

    if image_file is not None:
        simputfile += f" ImageFile={image_file.resolve()}"

    run_command(cmd=simputfile)


def simputmerge(infiles: list[Path], outfile: Path, fetch_extension: bool = True) -> None:
    ins = [str(f.resolve()) for f in infiles]
    ins = ",".join(ins)

    simputmerge = (
        f"simputmerge Infiles={ins} Outfile={outfile.resolve()} " f"FetchExtension={'yes' if fetch_extension else 'no'}"
    )

    run_command(cmd=simputmerge)
