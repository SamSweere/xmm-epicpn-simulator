import os
from pathlib import Path

import numpy as np
from astropy.io import fits

from src.xmm_utils.external_run import run_command


def compress_gzip(in_file_path: Path, out_file_path: Path, compresslevel=6, remove_file: bool = False):
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


def decompress_targz(in_file_path: Path, out_file_dir: Path):
    out_file_dir.mkdir(parents=True, exist_ok=True)
    run_command(f"tar -xzf {in_file_path.resolve()} -C {out_file_dir.resolve()} --strip-components=1")


def filter_event_pattern(eventlist_path: Path, max_event_pattern: int):
    if max_event_pattern == -1:
        # Use all event patterns
        return eventlist_path

    # Load the eventlist
    with fits.open(eventlist_path) as hdu:
        # Load the event header and data
        events_header = hdu["EVENTS"].header
        events_data = np.array(hdu["EVENTS"].data)

        # Filter the events data, remove every event with pattern type > max pattern type
        filtered_events_data = events_data[events_data["TYPE"] <= max_event_pattern]

        # Since we filtered the events, set the patterns to 0
        for i in range(max_event_pattern + 1, 13):
            events_header.set(f"NGRAD{i}", 0)
            events_header.set(f"NPGRA{i}", 0)

        events_header.set("NAXIS2", len(filtered_events_data))

        # Create new binary table from the filtered events
        filtered_events = fits.BinTableHDU(data=filtered_events_data, header=events_header)

        # Update the events with the new binary table
        hdu["EVENTS"] = filtered_events
        # Update history
        hdu["PRIMARY"].header["HISTORY"] = f"Removed all events with pattern type > {max_event_pattern}"

        # Overwrite the old eventlist
        hdu.writeto(eventlist_path, overwrite=True, checksum=True)

    return eventlist_path
