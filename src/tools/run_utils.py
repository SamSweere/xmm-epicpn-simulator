import os
from datetime import timedelta
from pathlib import Path

from loguru import logger

from src.config import XMM


def _opener(file, flags):
    fd = os.open(file, flags)
    os.chmod(fd, 0o777)
    return fd


def configure_logger(
    log_dir: Path,
    log_name: str,
    enqueue: bool,
    debug: bool,
    verbose: bool,
    retention: int,
    rotation: timedelta,
):
    if debug:
        log_level = "DEBUG"
    elif verbose:
        log_level = "INFO"
    else:
        log_level = "SUCCESS"
    log_file = log_dir / log_name
    logger.add(
        f"{log_file.resolve()}",
        enqueue=enqueue,
        level=log_level,
        rotation=rotation,
        retention=retention,
        opener=_opener,
    )


def load_satellites(instruments: dict) -> list:
    loaded_satellites = []

    if "xmm" in instruments:
        loaded_satellites.append(XMM(**instruments.pop("xmm")))

    if not loaded_satellites:
        raise RuntimeError()

    return loaded_satellites
