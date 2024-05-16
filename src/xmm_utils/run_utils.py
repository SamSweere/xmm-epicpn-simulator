import os
from datetime import timedelta
from pathlib import Path

from loguru import logger

from src.config import XMMInstrument


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


def load_satellites(instruments: dict) -> dict[dict]:
    loaded_satellites = {}

    if "xmm" in instruments:
        xmm = {}
        xmm_instruments: dict[str, dict] = XMMInstrument(**instruments.pop("xmm"))
        if xmm_instruments.emos1.use:
            xmm["emos1"] = xmm_instruments.emos1
        if xmm_instruments.emos2.use:
            xmm["emos2"] = xmm_instruments.emos2
        if xmm_instruments.epn.use:
            xmm["epn"] = xmm_instruments.epn

        if xmm:
            loaded_satellites["xmm"] = xmm

    if not loaded_satellites:
        raise RuntimeError()

    return loaded_satellites
