import os
import shutil
from pathlib import Path

import heasoftpy as hsp
from astropy.io import fits
from loguru import logger

from src.tools.external_run import run_command


def compress_gzip(in_file_path: Path, out_file_path: Path, compresslevel=6, remove_file: bool = False):
    out_file_path.parent.mkdir(exist_ok=True, parents=True)
    run_command(f"gzip -{compresslevel} -c {in_file_path.resolve()} > {out_file_path.resolve()}")
    if remove_file:
        in_file_path.unlink()


def compress_targz(in_path: Path, out_file_path: Path, remove_files: bool = False):
    if not out_file_path.name.endswith(".tar.gz"):
        raise ValueError(f"Output file path {out_file_path.resolve()} does not end with '.tar.gz'")
    out_file_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = " --remove-files" if remove_files else ""
    run_command(
        f"cd {in_path.parent.resolve()} && "
        + f"tar -czf {out_file_path.resolve()} {in_path.name}{os.sep} --overwrite{suffix}"
    )
    if remove_files:
        shutil.rmtree(in_path)


def decompress_targz(in_file_path: Path, out_file_dir: Path, tar_options: str = ""):
    out_file_dir.mkdir(parents=True, exist_ok=True)
    run_command(f"tar -xzf {in_file_path.resolve()} -C {out_file_dir.resolve()} {tar_options}")
    logger.success(f"Decompressed {in_file_path} to {out_file_dir}")


def filter_event_pattern(eventlist_path: Path, max_event_pattern: int) -> Path | None:
    if max_event_pattern == -1 or max_event_pattern == 12:
        # Use all event patterns
        logger.debug(f"There is nothing to filter for {eventlist_path}.")
        return eventlist_path

    logger.debug(f"Filtering {eventlist_path} for pattern <= {max_event_pattern}.")

    outfile = eventlist_path.parent / f"{eventlist_path.stem}_filtered.fits"

    with hsp.utils.local_pfiles_context():
        # Filter events
        hsp.ftcopy(
            infile=f"{eventlist_path}[EVENTS][TYPE <= {max_event_pattern}]",
            outfile=f"{outfile}",
            history="yes",
        )

        assert outfile.exists()

        eventlist_path.unlink()

        infile = f"{outfile}[EVENTS]"
        for i in range(max_event_pattern + 1, 13):
            hsp.fthedit(infile=infile, keyword=f"NGRAD{i}", operation="add", value=0)
            hsp.fthedit(infile=infile, keyword=f"NPGRA{i}", operation="add", value=0)

        data = fits.getdata(outfile, "EVENTS")
        size = data.size
        del data

        if size == 0:
            # No events left after filtering
            outfile.unlink()
            return None

    return outfile
