import subprocess

from loguru import logger


def run_command(cmd, cmd_input: str | None = None) -> None:
    logger.debug(f"Running command:\n\t{cmd}")
    # Execute a shell command with the stdout and stderr being redirected to a log file
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
