import os
from datetime import timedelta
from pathlib import Path
from typing import List

from loguru import logger


def opener(file, flags):
    return os.open(file, flags, 0o777)


def configure_logger(
        log_dir: Path,
        log_name: str,
        enqueue: bool,
        debug: bool,
        retention: int,
        rotation: timedelta
):
    log_level = "DEBUG" if debug else "INFO"
    log_file = log_dir / log_name
    logger.add(log_file.resolve(), enqueue=enqueue, level=log_level, rotation=rotation, retention=retention,
               opener=opener)


def handle_error(error):
    logger.exception(error)


def create_dirs(
        list_of_paths: List[Path]
) -> None:
    for path in list_of_paths:
        path.mkdir(parents=True, exist_ok=True)
