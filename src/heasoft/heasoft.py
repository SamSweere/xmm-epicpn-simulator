import os
import subprocess
from os.path import exists, isfile, join
from pathlib import Path
from tempfile import TemporaryDirectory

from loguru import logger


def _run_cmd(cmd: list[str], timeout: float = 3600) -> None:
    with TemporaryDirectory(
        prefix="heasoft_",
    ) as tmpdir:
        cmd.extend(["history=yes", "chatter=1"])
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


def ftmerge(
    infiles: list[Path],
    outfile: Path,
    consume_data: bool,
) -> Path:
    exec_cmd = join(os.environ["HEADAS"], "bin", "ftmerge")

    assert exists(exec_cmd) and isfile(exec_cmd), f"{exec_cmd} does not exist"

    for infile in infiles:
        assert exists(infile)

    cmd_params = [f"infile={','.join([str(infile) for infile in infiles])}", f"outfile={outfile}", "clobber=yes"]
    cmd = [exec_cmd, *cmd_params]
    _run_cmd(cmd)

    if consume_data:
        for infile in infiles:
            infile.unlink()

    assert exists(outfile)

    return outfile


def ftcopy(
    infile: str | Path,
    outfile: Path,
) -> Path:
    exec_cmd = join(os.environ["HEADAS"], "bin", "ftcopy")

    assert exists(exec_cmd)
    if isinstance(infile, Path):
        assert exists(infile)

    if isinstance(infile, str):
        assert exists(infile.split("[")[0])

    cmd_params = [f"infile={infile}", f"outfile={outfile}", "clobber=yes"]
    cmd = [exec_cmd, *cmd_params]
    _run_cmd(cmd)

    assert exists(outfile)

    return outfile


def fthedit(
    infile: str | Path,
    keyword: str,
    operation: str,
    value: str,
    **kwargs,
) -> None:
    exec_cmd = join(os.environ["HEADAS"], "bin", "fthedit")

    assert exists(exec_cmd)
    if isinstance(infile, Path):
        assert exists(infile)

    if isinstance(infile, str):
        assert exists(infile.split("[")[0])

    cmd_params = [f"infile={infile}", f"keyword={keyword}", f"operation={operation}", f"value={value}"]
    for k, v in kwargs.items():
        cmd_params.append(f"{k}={v}")
    cmd = [exec_cmd, *cmd_params]
    _run_cmd(cmd)


def ftedit(
    infile: str | Path,
    column: str,
    row: int,
    value: str,
) -> None:
    exec_cmd = join(os.environ["HEADAS"], "bin", "ftedit")

    assert exists(exec_cmd)
    if isinstance(infile, Path):
        assert exists(infile)

    if isinstance(infile, str):
        assert exists(infile.split("[")[0])

    cmd_params = [f"infile={infile}", f"column={column}", f"row={row}", f"value={value}"]
    cmd = [exec_cmd, *cmd_params]
    _run_cmd(cmd)


def ftimgcalc(
    outfile: Path,
    expr: str,
    **kwargs,
) -> Path:
    exec_cmd = join(os.environ["HEADAS"], "bin", "ftimgcalc")

    assert exists(exec_cmd)

    cmd_params = [f"outfile={outfile}", f"expr={expr}", "clobber=yes"]
    for k, v in kwargs.items():
        cmd_params.append(f"{k}={v}")

    cmd = [exec_cmd, *cmd_params]
    _run_cmd(cmd)

    assert exists(outfile)

    return outfile
