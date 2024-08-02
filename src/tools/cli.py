import subprocess

from loguru import logger


def run_command(cmd, cmd_input: str | None = None) -> None:
    logger.debug(f"Running command:\n\t{cmd}")
    # Execute a shell command with the stdout and stderr being redirected to a log file
    cmd = [cmd]
    if cmd_input is not None:
        cmd_input = bytes(cmd_input, "ascii")
    result = subprocess.run(
        args=cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        input=cmd_input,
        close_fds=True,
    )
    retcode = result.returncode
    if retcode != 0:
        raise RuntimeError(f"Execution of {cmd} returned {retcode}\n{result.stdout}")


def run(cmd: list[str], timeout: float = 3600) -> None:
    proc = subprocess.run(
        cmd,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    try:
        returncode = proc.returncode

        if returncode != 0 or proc.stderr:
            logger.error(
                f"{' '.join(cmd)} failed with return code '{returncode}'\n"
                f"\tand output '{proc.stdout}'\n"
                f"\tand stderr '{proc.stderr}'"
            )
            raise subprocess.CalledProcessError(returncode, cmd, output=proc.stdout, stderr=proc.stderr)

        logger.debug(f"{' '.join(cmd)} output: {proc.stdout}")
    except subprocess.TimeoutExpired as e:
        logger.error(
            f"{' '.join(cmd)} timed out after '{e.timeout}'\n"
            f"\tand output '{e.output}'\n"
            f"\tand stderr '{e.stderr}'"
        )
        raise
