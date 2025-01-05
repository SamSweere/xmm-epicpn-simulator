from pathlib import Path

from astropy.io import fits
from loguru import logger

import src.heasoft as hsp


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
                clobber="yes",
            )

            assert outfile.exists()

            for ext in ["PRIMARY", "EVENTS", "STDGTI"]:
                hsp.fthedit(
                    infile=f"{outfile}[{ext}]",
                    keyword="TSTART",
                    operation="add",
                    value=f"{t_start}",
                    unit="s",
                )
                hsp.fthedit(
                    infile=f"{outfile}[{ext}]",
                    keyword="TSTOP",
                    operation="add",
                    value=f"{t_stop}",
                    unit="s",
                )

            hsp.fthedit(
                infile=f"{outfile}[EVENTS]",
                keyword="EXPOSURE",
                operation="add",
                value=f"{split}",
                unit="s",
            )
            hsp.ftedit(
                infile=f"{outfile}[STDGTI]",
                column="START",
                row=1,
                value=f"{t_start}",
            )
            hsp.ftedit(
                infile=f"{outfile}[STDGTI]",
                column="STOP",
                row=1,
                value=f"{t_stop}",
            )

            split_exps.append((outfile, split))

    if consume_data:
        eventlist_path.unlink()

    return split_exps
