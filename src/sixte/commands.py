import os
import subprocess
from os.path import exists, isfile, join
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal

from loguru import logger


def _run_cmd(cmd: list[str], timeout: float = 3600) -> None:
    with TemporaryDirectory(prefix="sixte_") as tmpdir:
        env = os.environ.copy()
        env["PFILES"] = f"{tmpdir}:{os.environ['PFILES']}"
        proc = subprocess.run(
            args=cmd,
            env=env,
            capture_output=True,
            text=True,
            timeout=timeout,
        )

    try:
        returncode = proc.returncode

        if returncode != 0 or proc.stderr:
            logger.error(
                f"{cmd} failed with return code '{returncode}'\n"
                f"\tand output '{proc.stdout}'\n"
                f"\tand stderr '{proc.stderr}'"
            )
            raise subprocess.CalledProcessError(returncode, cmd, output=proc.stdout, stderr=proc.stderr)

        logger.debug(f"{cmd} output: {proc.stdout}")
    except subprocess.TimeoutExpired as e:
        logger.error(
            f"{cmd} timed out after '{e.timeout}'\n" f"\tand output '{e.output}'\n" f"\tand stderr '{e.stderr}'"
        )
        raise


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
    chatter: int = 1,
    history: bool = True,
) -> Path:
    exec_cmd = join(os.environ["SIXTE"], "bin", "imgev")

    assert exists(exec_cmd) and isfile(exec_cmd), f"{exec_cmd} does not exist"
    assert exists(evt_file) and isfile(evt_file), f"{evt_file} does not exist"

    cmd_params = [
        f"EvtFile={evt_file}",
        f"Image={image}",
        f"CoordinateSystem={coordinate_system}",
        f"NAXIS1={naxis1}",
        f"NAXIS2={naxis2}",
        f"CUNIT1={cunit1}",
        f"CUNIT2={cunit2}",
        f"CRVAL1={crval1}",
        f"CRVAL2={crval2}",
        f"CRPIX1={crpix1}",
        f"CRPIX2={crpix2}",
        f"CDELT1={cdelt1}",
        f"CDELT2={cdelt2}",
        f"chatter={chatter}",
        f"clobber={str(clobber).lower()}",
        f"history={str(history).lower()}",
    ]

    cmd = [exec_cmd, *cmd_params]
    _run_cmd(cmd)

    assert exists(image) and isfile(image), f"{image} does not exist"

    return image


def sixtesim(
    output_path: Path,
    xml_file: Path,
    ra: float,
    dec: float,
    rollangle: float,
    simput: Path,
    exposure: int,
    chatter: int = 1,
    clobber: bool = True,
    history: bool = True,
) -> None:
    exec_cmd = join(os.environ["SIXTE"], "bin", "sixtesim")

    assert exists(exec_cmd) and isfile(exec_cmd), f"{exec_cmd} does not exist"
    assert exists(xml_file) and isfile(xml_file), f"{xml_file} does not exist"
    assert exists(simput) and isfile(simput), f"{simput} does not exist"

    cmd_params = [
        f"Prefix={output_path}/",
        f"XMLFile={xml_file}",
        f"RA={ra}",
        f"DEC={dec}",
        f"rollangle={rollangle}",
        f"Simput={simput}",
        f"Exposure={exposure}",
        f"chatter={chatter}",
        f"clobber={str(clobber).lower()}",
        "Background=no",
        f"history={str(history).lower()}",
        "progressbar=n",
    ]

    cmd = [exec_cmd, *cmd_params]
    _run_cmd(cmd, 1.1 * exposure)

    with os.scandir(output_path) as it:
        assert any(it)


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
    exec_cmd = join(os.environ["SIXTE"], "bin", "simputfile")

    assert exists(exec_cmd) and isfile(exec_cmd), f"{exec_cmd} does not exist"

    cmd_params = [
        f"Simput={simput}",
        f"RA={ra}",
        f"DEC={dec}",
        f"srcFlux={src_flux}",
        f"Emin={emin}",
        f"Emax={emax}",
        f"clobber={str(clobber).lower()}",
        f"chatter={chatter}",
        f"history={str(history).lower()}",
    ]

    if xspec_file is not None:
        assert exists(xspec_file) and isfile(xspec_file), f"{xspec_file} does not exist"
        cmd_params.append(f"XSPECFile={xspec_file}")

    if ascii_file is not None:
        assert exists(ascii_file) and isfile(ascii_file), f"{ascii_file} does not exist"
        cmd_params.append(f"ASCIIFile={ascii_file}")

    if image_file is not None:
        assert exists(image_file) and isfile(image_file), f"{image_file} does not exist"
        cmd_params.append(f"ImageFile={image_file}")

    cmd = [exec_cmd, *cmd_params]
    _run_cmd(cmd)


def simputmerge(infiles: list[Path], outfile: Path, fetch_extension: bool = True) -> None:
    exec_cmd = join(os.environ["SIXTE"], "bin", "simputmerge")

    assert exists(exec_cmd) and isfile(exec_cmd), f"{exec_cmd} does not exist"

    for infile in infiles:
        assert exists(infile) and isfile(infile), f"{infile} does not exist"

    cmd_params = [
        f"Infiles={','.join([str(infile) for infile in infiles])}",
        f"Outfile={outfile.resolve()}",
        f"FetchExtension={'yes' if fetch_extension else 'no'}",
    ]

    cmd = [exec_cmd, *cmd_params]
    _run_cmd(cmd)
