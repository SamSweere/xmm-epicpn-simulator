import subprocess
from typing import Optional

from loguru import logger


def _run_command(cmd, verbose=True, cmd_input: Optional[str] = None) -> None:
    if verbose:
        logger.info(f"Running command:\n\t{cmd}")
    #
    # Execute a shell command with the stdout and stderr being redirected to a log file
    #
    if cmd_input is not None:
        cmd_input = bytes(cmd_input, 'ascii')
    result = subprocess.run(args=cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, input=cmd_input,
                            close_fds=True)
    retcode = result.returncode
    if retcode != 0:
        raise RuntimeError(f"Execution of {cmd} returned {retcode}\n{result.stdout}")


def run_headas_command(
        cmd: str,
        cmd_input: Optional[str] = None,
        verbose: bool = True
):
    _run_command(cmd, cmd_input=cmd_input, verbose=verbose)
