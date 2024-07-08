from pathlib import Path

import heasoftpy as hsp
from astropy.io import fits
from loguru import logger


def merge_ccd_eventlists(infiles: list[Path], out_dir: Path, consume_data: bool) -> Path:
    if len(infiles) == 1:
        return infiles[0]

    # See https://www.sternwarte.uni-erlangen.de/research/sixte/data/simulator_manual_v1.3.11.pdf for information
    all_ccds = [f"{infile.resolve()}" for infile in infiles]
    outfile = out_dir / "merged.fits"
    params = {
        "infile": ",".join(all_ccds),
        "outfile": f"{outfile.resolve()}",
        "clobber": "yes",
    }

    with hsp.utils.local_pfiles_context():
        hsp.ftmerge(params)

    logger.success(f"Successfully ran 'ftmerge' with params: {params}")

    if consume_data:
        for infile in infiles:
            infile.unlink()

    return outfile


def split_eventlist(run_dir: Path, eventlist_path: Path, consume_data: bool, multiples: int = 10000) -> list[Path]:
    # This function splits an eventlist in multiples of multiples and saves them.
    # It returns the split files
    logger.debug(f"Splitting {eventlist_path}")
    exposure = int(fits.getheader(eventlist_path, "EVENTS")["EXPOSURE"])
    split_exps = []

    with hsp.utils.local_pfiles_context():
        for split in range(multiples, exposure + multiples, multiples):
            num = int(exposure / split)
            for i in range(num):
                t_start = i * split
                t_stop = (i + 1) * split

                outfile = run_dir / f"{round(split / 1000)}ks_p_{i}-{num - 1}_evt.fits"

                hsp.ftcopy(
                    infile=f"{eventlist_path}[EVENTS][TIME >= {t_start} && TIME < {t_stop}]",
                    outfile=f"{outfile}",
                    clobber="yes",
                    copyall="yes",
                )

                assert outfile.exists()

                hsp.fthedit(
                    infile=f"{outfile}[PRIMARY]",
                    keyword="TSTART",
                    operation="add",
                    value=t_start,
                )
                hsp.fthedit(
                    infile=f"{outfile}[PRIMARY]",
                    keyword="TSTOP",
                    operation="add",
                    value=t_stop,
                )
                hsp.fthedit(
                    infile=f"{outfile}[EVENTS]",
                    keyword="TSTART",
                    operation="add",
                    value=t_start,
                )
                hsp.fthedit(
                    infile=f"{outfile}[EVENTS]",
                    keyword="TSTOP",
                    operation="add",
                    value=t_stop,
                )
                hsp.fthedit(
                    infile=f"{outfile}[EVENTS]",
                    keyword="EXPOSURE",
                    operation="add",
                    value=split,
                )
                hsp.fthedit(
                    infile=f"{outfile}[STDGTI]",
                    keyword="TSTART",
                    operation="add",
                    value=t_start,
                )
                hsp.fthedit(
                    infile=f"{outfile}[STDGTI]",
                    keyword="TSTOP",
                    operation="add",
                    value=t_stop,
                )
                hsp.ftedit(
                    infile=f"{outfile}[STDGTI]",
                    column="START",
                    row=1,
                    value=t_start,
                )
                hsp.ftedit(
                    infile=f"{outfile}[STDGTI]",
                    column="STOP",
                    row=1,
                    value=t_stop,
                )
                split_exps.append(outfile)

    if consume_data:
        eventlist_path.unlink()

    return split_exps
