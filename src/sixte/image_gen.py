from pathlib import Path
from typing import List

import heasoftpy as hsp
from astropy.io import fits

from src.xmm_utils.external_run import run_headas_command


def merge_ccd_eventlists(
        infiles: List[Path],
        out_dir: Path,
        keep_files: bool = False,
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
    with hsp.utils.local_pfiles_context():
        hsp.ftmerge(params)

    if verbose:
        print(f"Successfully ran 'ftmerge' with params: {params}")

    if not keep_files:
        for ccd in all_ccds:
            Path(ccd).unlink()

    return outfile


# Merge and generate image command
def imgev(
        evt_file: Path,
        input_folder: Path,
        out_name: str,
        naxis1: int,
        naxis2: int,
        crval1: float,
        crval2: float,
        crpix1: float,
        crpix2: float,
        cdelt1: float,
        cdelt2: float,
        keep_files: bool = False,
        verbose: bool = True
) -> Path:
    outfile = input_folder / out_name
    cmd = (f"imgev EvtFile='{evt_file.resolve()}' Image='{outfile.resolve()}' CoordinateSystem=0 "
           f"Projection=TAN NAXIS1='{naxis1}' NAXIS2='{naxis2}' CUNIT1=deg CUNIT2=deg CRVAL1='{crval1}' "
           f"CRVAL2='{crval2}' CRPIX1='{crpix1}' CRPIX2='{crpix2}' CDELT1='{cdelt1}' CDELT2='{cdelt2}' clobber=yes")

    run_headas_command(cmd=cmd, verbose=verbose)

    if not keep_files:
        evt_file.unlink()

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
