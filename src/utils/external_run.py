import subprocess
from pathlib import Path
from typing import List, Optional

from utils import log


def _run_command(command, verbose=True, input=None):
    if verbose:
        print(command)
    #
    # Execute a shell command with the stdout and stderr being redirected to a log file
    #
    # shell=False suggested Francesco Pierfederici <fra.pierfederici@icloud.com>, but it does not work as expected
    #
    # not using the other suggestion check=True as I am not sure it will do what I need. It will raise an exception that I am
    # catching with try: except: anyway.
    #
    retcode = None
    try:
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, input=input,
                                encoding='ascii')
        retcode = result.returncode
        if retcode < 0:
            w_message = f"Execution of {command} was terminated by signal: {-retcode} \n {result.stdout}"
            log.wlog(w_message, verbose=verbose)
        elif retcode > 0:
            message = f"Execution of {command} returned {retcode}, \n {retcode, result.stdout}"
            log.elog(message, verbose=verbose)
            raise RuntimeError(message)
        else:
            message = f"Execution of {command} returned {retcode}, \n {retcode, result.stdout}"
            log.plog(message, verbose=verbose)
    except OSError as e:
        e_message = f"Execution of {command} failed: {e}"
        log.elog(e_message)
        raise OSError(e_message)
    return retcode, result


def run_headas_command(
        command: str,
        run_dir: Optional[Path] = None,
        input=None,
        verbose=True
):
    # Initilialize HEADAS before running the command
    cmds: List[str] = ["" if run_dir is None else f"cd {run_dir.resolve()}",
                       ". $HEADAS/headas-init.sh",
                       "export HEADASNOQUERY= && export HEADASPROMPT=/dev/null",
                       command]
    cmd = " && ".join(cmds)
    _run_command(cmd, input=input, verbose=verbose)
