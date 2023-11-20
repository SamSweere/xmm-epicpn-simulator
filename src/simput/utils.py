import os
import random
import shutil
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
from typing import List, Dict, Literal, Optional

from loguru import logger
from xspec import Model, Xset

from src.sixte import commands


def get_spectrumfile(run_dir: Path, norm=0.01, verbose=True):
    spectrum_file = run_dir / "spectrum.xcm"
    if not spectrum_file.exists():
        if verbose:
            logger.info(f"Spectrum file at {spectrum_file.resolve()} does not exist. Will create a new one.")

        with open(os.devnull, "w") as devnull, redirect_stdout(devnull), redirect_stderr(devnull):
            Model("phabs*power", setPars={1: 0.04, 2: 2.0, 3: norm})
            Xset.save(f"{spectrum_file.resolve()}")

    return spectrum_file


def merge_simputs(
        simput_files: List[Path],
        output_file: Path,
        verbose=True
) -> Path:
    # Combine the simput point sources
    if len(simput_files) == 1:
        file = simput_files[0]
        shutil.copy2(file, output_file)
    else:
        commands.simputmerge(infiles=simput_files, outfile=output_file, fetch_extension=True)

    return output_file


def _order(iterable: List, order: str) -> List:
    if order == "normal":
        pass
    elif order == "reversed":
        iterable.reverse()
    elif order == "random":
        random.shuffle(iterable)
    else:
        raise ValueError(f'Order: {order} not in known options list of "normal", "reversed", "random"')

    return iterable


def get_simputs(
        instrument_name: Optional[Literal["epn", "emos1", "emos2"]],
        simput_path: Path,
        mode_amount_dict: Dict[str, int],
        order='normal'
) -> Dict[str, List[Path]]:
    # Order options: normal (front to back), reversed (back to front), random

    simput_files: Dict[str, List[Path]] = {}
    for mode, amount in mode_amount_dict.items():
        if amount == 0:
            continue

        if mode == "background":
            mode_path = simput_path / instrument_name / mode
        else:
            mode_path = simput_path / mode
        files = mode_path.glob("*.simput.gz")

        if amount == -1:
            # Do all
            simput_files[mode] = _order(list(files), order)
        else:
            tmp = []
            for file_count, file in enumerate(files):
                if file_count < amount:
                    tmp.append(file)
                else:
                    break
            simput_files[mode] = _order(tmp, order)

    return simput_files
