import os
import shutil
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from loguru import logger
from xspec import Model, Xset

from src.sixte import commands


def get_spectrumfile(run_dir: Path, norm=0.01) -> Path:
    spectrum_file = run_dir / "spectrum.xcm"
    if not spectrum_file.exists():
        logger.info(f"Spectrum file at {spectrum_file.resolve()} does not exist. Will create a new one.")

        with open(os.devnull, "w") as f, redirect_stdout(f), redirect_stderr(f):
            Model("phabs*power", setPars={1: 0.04, 2: 2.0, 3: norm})
            Xset.save(f"{spectrum_file.resolve()}")

    return spectrum_file


def merge_simputs(simput_files: list[Path], output_file: Path) -> Path:
    # Combine the simput point sources
    if len(simput_files) == 1:
        file = simput_files[0]
        shutil.copy2(file, output_file)
    else:
        commands.simputmerge(infiles=simput_files, outfile=output_file, fetch_extension=True)

    return output_file
