from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List

import heasoftpy as hsp
from astropy.io import fits
from loguru import logger


def merge_ccd_eventlists(
        infiles: List[Path],
        out_dir: Path,
        verbose: bool = True
) -> Path:
    # See https://www.sternwarte.uni-erlangen.de/research/sixte/data/simulator_manual_v1.3.11.pdf for information
    all_ccds = [f"{infile.resolve()}" for infile in infiles]
    outfile = out_dir / "merged.fits"
    params = {
        "infile": ",".join(all_ccds),
        "outfile": f"{outfile.resolve()}",
        "clobber": "yes"
    }

    with TemporaryDirectory(prefix="hsp_") as tmp_dir, hsp.utils.local_pfiles_context(tmp_dir):
        hsp.ftmerge(params)

    if verbose:
        logger.info(f"Successfully ran 'ftmerge' with params: {params}")

    return outfile


def split_eventlist(run_dir: Path, eventlist_path: Path, multiples: int = 10000, verbose: bool = True):
    # This function splits an eventlist in multiples of multiples and saves them.
    # It returns the split files
    with fits.open(eventlist_path, mode="readonly") as hdu:
        exposure = int(hdu['EVENTS'].header['EXPOSURE'])
        split_exposure_evt_files = []

        events_data = hdu['EVENTS'].data.copy()

        for split_exp in range(multiples, exposure + multiples, multiples):
            num = int(exposure / split_exp)
            # print(f"{num} x split exposure {split_exp} s")

            for i in range(num):
                t_start = i * split_exp
                t_stop = (i + 1) * split_exp

                # Filter the data
                mask = events_data["TIME"] >= t_start
                mask = (mask == (events_data["TIME"] < t_stop))
                hdu["EVENTS"].data = events_data[mask]

                # Update the header
                hdu['PRIMARY'].header['TSTART'] = t_start
                hdu['PRIMARY'].header['TSTOP'] = t_stop

                hdu['EVENTS'].header['TSTART'] = t_start
                hdu['EVENTS'].header['TSTOP'] = t_stop
                hdu['EVENTS'].header['EXPOSURE'] = split_exp

                hdu['STDGTI'].header['TSTART'] = t_start
                hdu['STDGTI'].header['TSTOP'] = t_stop

                hdu['STDGTI'].data[0] = (float(t_start), float(t_stop))

                base_name = f"{round(split_exp / 1000)}ks_p_{i}-{num - 1}"
                outfile = run_dir / f"{base_name}_evt.fits"
                split_exposure_evt_files.append({'outfile': outfile.resolve(),
                                                 'base_name': base_name,
                                                 't_start': t_start,
                                                 't_stop': t_stop,
                                                 'split_num': i,
                                                 'total_splits': num,
                                                 'exposure': split_exp})

                hdu.writeto(outfile)

    return split_exposure_evt_files
