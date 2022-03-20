# Generate the final fits image
import os

from astropy.io import fits

from utils.external_run import run_headas_command


def run_merge_command(input_folder, command):
    # Run the command in the input folder
    merge_command = f"cd {input_folder} && {command}"

    run_headas_command(merge_command)


def merge_ccd_eventlists(input_folder):
    final_evt_name = 'ccd_combined_evt.fits'

    commands = []
    # Run the command in the input folder
    # commands.append(f"cd {input_folder}")

    # See https://www.sternwarte.uni-erlangen.de/research/sixte/data/simulator_manual_v1.3.11.pdf for information
    # Merge the ccds of every quadrant
    commands.append('ftmerge "ccd01_evt.fits","ccd02_evt.fits","ccd03_evt.fits" "ccd_q0_evt.fits" clobber=yes')
    commands.append('ftmerge "ccd04_evt.fits","ccd05_evt.fits","ccd06_evt.fits" "ccd_q1_evt.fits" clobber=yes')
    commands.append('ftmerge "ccd07_evt.fits","ccd08_evt.fits","ccd09_evt.fits" "ccd_q2_evt.fits" clobber=yes')
    commands.append('ftmerge "ccd10_evt.fits","ccd11_evt.fits","ccd12_evt.fits" "ccd_q3_evt.fits" clobber=yes')

    # Merge the quadrants into one
    commands.append(f'ftmerge "ccd_q0_evt.fits","ccd_q1_evt.fits","ccd_q2_evt.fits","ccd_q3_evt.fits" '
                    f'"{final_evt_name}" clobber=yes')  # chatter=5

    for command in commands:
        run_merge_command(input_folder=input_folder, command=command)

    # Return the final event path
    return os.path.join(input_folder, final_evt_name)


# %%
# Merge and generate image command
def generate_fits_image(evt_file, input_folder, out_name, naxis1, naxis2, crval1, crval2, crpix1, crpix2, cdelt1,
                        cdelt2, verbose=True):
    commands = []
    # Run the command in the input folder
    commands.append(f"cd {input_folder}")

    # Generate the image
    commands.append(f"imgev EvtFile='{evt_file}' Image='{out_name}' CoordinateSystem=0 Projection=TAN "
                    f"NAXIS1='{naxis1}' NAXIS2='{naxis2}' CUNIT1=deg CUNIT2=deg CRVAL1='{crval1}' CRVAL2='{crval2}' "
                    f"CRPIX1='{crpix1}' CRPIX2='{crpix2}' CDELT1='{cdelt1}' CDELT2='{cdelt2}' clobber=yes")

    # Combine the commands into one string
    merge_gen_command = " && ".join(commands)

    # merge_gen_command = f"/home/sam/Documents/ESA/xmm-superres/simulation/ftmerge_imgev.sh " \
    #                     f"'{input_folder}' '{out_name}' '{naxis1}' '{naxis2}' '{crval1}' '{crval2}' '{crpix1}' '{crpix2}' " \
    #                     f" '{cdelt1}' '{cdelt2}' '{dtype}'"
    run_headas_command(merge_gen_command, verbose=verbose)


def split_eventlist(run_dir, eventlist_path, multiples=10000, verbose=True):
    # This function splits an eventlist in multiples of multiples and saves them.
    # It returns the split files
    hdu = fits.open(eventlist_path)

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

            hdu.writeto(os.path.join(run_dir, filename))

    return split_exposure_evt_files
