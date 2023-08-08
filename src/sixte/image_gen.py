# Generate the final fits image
from pathlib import Path

from astropy.io import fits
from tqdm import tqdm

from src.utils.external_run import run_headas_command


def merge_ccd_eventlists(input_folder: Path, out_filename: str, verbose: bool = False):
    # Run the command in the input folder
    # commands.append(f"cd {input_folder}")

    # See https://www.sternwarte.uni-erlangen.de/research/sixte/data/simulator_manual_v1.3.11.pdf for information
    # Merge the ccds of every quadrant
    # TODO fix path to individual .fits files
    cmds = ['ftmerge "ccd01_evt.fits","ccd02_evt.fits","ccd03_evt.fits" "ccd_q0_evt.fits" clobber=yes',
            'ftmerge "ccd04_evt.fits","ccd05_evt.fits","ccd06_evt.fits" "ccd_q1_evt.fits" clobber=yes',
            'ftmerge "ccd07_evt.fits","ccd08_evt.fits","ccd09_evt.fits" "ccd_q2_evt.fits" clobber=yes',
            'ftmerge "ccd10_evt.fits","ccd11_evt.fits","ccd12_evt.fits" "ccd_q3_evt.fits" clobber=yes',
            f'ftmerge "ccd_q0_evt.fits","ccd_q1_evt.fits","ccd_q2_evt.fits","ccd_q3_evt.fits" '
            f'"{out_filename}" clobber=yes']

    # Merge the quadrants into one
    for cmd in tqdm(cmds, desc="Running commands to merge CCD eventlists"):
        run_headas_command(cmd=cmd, verbose=verbose)

    # Return the final event path
    return input_folder / out_filename


# Merge and generate image command
def generate_fits_image(evt_file, input_folder, out_name, naxis1, naxis2, crval1, crval2, crpix1, crpix2, cdelt1,
                        cdelt2, verbose=True):
    cmds = [
        f"imgev EvtFile='{input_folder / evt_file}' Image='{input_folder / out_name}' CoordinateSystem=0 Projection=TAN "
        f"NAXIS1='{naxis1}' NAXIS2='{naxis2}' CUNIT1=deg CUNIT2=deg CRVAL1='{crval1}' CRVAL2='{crval2}' "
        f"CRPIX1='{crpix1}' CRPIX2='{crpix2}' CDELT1='{cdelt1}' CDELT2='{cdelt2}' clobber=yes"]

    for cmd in tqdm(cmds, desc="Generating FITS files"):
        run_headas_command(cmd=cmd, verbose=verbose)


def split_eventlist(run_dir: Path, eventlist_path: Path, multiples: int = 10000, verbose: bool = True):
    # This function splits an eventlist in multiples of multiples and saves them.
    # It returns the split files
    with fits.open(eventlist_path, mode="readonly") as hdu:
        exposure = int(hdu['EVENTS'].header['EXPOSURE'])
        split_exposure_evt_files = []

        for split_exp in range(multiples, exposure + multiples, multiples):
            num = int(exposure / split_exp)
            # print(f"{num} x split exposure {split_exp} s")

            for i in range(num):
                t_start = i * split_exp
                t_stop = (i + 1) * split_exp
                # print(f"Time range: ({t_start},{t_stop})")

                # Load the full eventlist file
                # TODO I think that this file opening is not needed
                hdu = fits.open(eventlist_path)

                # Filter the data
                data = hdu['EVENTS'].data
                mask = data['TIME'] >= t_start
                mask = (mask == (data['TIME'] < t_stop))
                hdu['EVENTS'].data = data[mask]

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
                filename = base_name + "_evt.fits"
                split_exposure_evt_files.append({'filename': filename,
                                                 'base_name': base_name,
                                                 't_start': t_start,
                                                 't_stop': t_stop,
                                                 'split_num': i,
                                                 'total_splits': num,
                                                 'exposure': split_exp})

                hdu.writeto(run_dir / filename)

    return split_exposure_evt_files
