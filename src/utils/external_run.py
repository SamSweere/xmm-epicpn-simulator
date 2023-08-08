import subprocess
from typing import Optional


def _run_command(cmd, verbose=True, cmd_input: Optional[str] = None) -> None:
    if verbose:
        print(f"Running command: {cmd}")
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
    # Initilialize HEADAS before running the command
    cmd = f". $HEADAS/headas-init.sh && {cmd}"
    _run_command(cmd, cmd_input=cmd_input, verbose=verbose)


def run_sixte_command(
        cmd: str,
        verbose: bool = True
) -> None:
    cmd = (f"export SIMPUT=/home/bojantodorkov/simput &&"
           f"export SIXTE=/home/bojantodorkov/simput &&"
           f". $SIXTE/bin/sixte-install.sh && {cmd}")
    run_headas_command(cmd, verbose=verbose)
