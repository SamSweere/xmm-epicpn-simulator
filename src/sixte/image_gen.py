from pathlib import Path

from astropy.io import fits
from loguru import logger

from src.heasoft import heasoft as hsp


def split_eventlist(
    run_dir: Path, eventlist_path: Path, consume_data: bool, multiples: int = 10000
) -> list[tuple[Path, int]]:
    # This function splits an eventlist in multiples of multiples and saves them.
    # It returns the split files
    logger.debug(f"Splitting {eventlist_path}")
    exposure = int(fits.getheader(eventlist_path, "EVENTS")["EXPOSURE"])
    split_exps = []

    for split in range(multiples, exposure + multiples, multiples):
        num = int(exposure / split)
        for i in range(num):
            t_start = i * split
            t_stop = (i + 1) * split

            outfile = run_dir / f"{round(split / 1000)}ks_p_{i}-{num - 1}_evt.fits"

            outfile = hsp.ftcopy(
                infile=f"{eventlist_path}[EVENTS][TIME >= {t_start} && TIME < {t_stop}]",
                outfile=outfile,
            )

            assert outfile.exists()

            for ext in ["PRIMARY", "EVENTS", "STDGTI"]:
                hsp.fthedit(f"{outfile}[{ext}]", "TSTART", "add", f"{t_start}", unit="s")
                hsp.fthedit(f"{outfile}[{ext}]", "TSTOP", "add", f"{t_stop}", unit="s")

            hsp.fthedit(f"{outfile}[EVENTS]", "EXPOSURE", "add", f"{split}", unit="s")
            hsp.ftedit(f"{outfile}[STDGTI]", "START", 1, f"{t_start}")
            hsp.ftedit(f"{outfile}[STDGTI]", "STOP", 1, f"{t_stop}")

            split_exps.append((outfile, split))

    if consume_data:
        eventlist_path.unlink()

    return split_exps
