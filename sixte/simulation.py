import os

from utils.external_run import run_headas_command


def run_sixte(name, simput_file_path, xml_file_path, output_dir, exposure, ra, dec, rollangle, verbose):
    if verbose:
        print("--------------------")
        print("Running " + name)

    raw_filepath = os.path.join(output_dir, name + "_raw.fits")
    evt_filepath = os.path.join(output_dir, name + "_evt.fits")

    sixte_command = f"runsixt RawData={raw_filepath} EvtFile={evt_filepath} Mission='XMM' Instrument='EPICPN' " \
                    f"Mode='FFTHIN' XMLFile={xml_file_path} Simput={simput_file_path} Exposure={exposure} RA={ra} " \
                    f"Dec={dec} rollangle={rollangle} clobber=yes"

    # Run the command as headas command (this will initialize headas before running this
    run_headas_command(sixte_command, verbose=verbose)


def run_sixte_all_ccds(xmm, simput_file_path, rundir, exposure, ra=0.0, dec=0.0, rollangle=0.0, verbose=True):
    # Run the sixte in paralell
    # pool = get_multiprocessing_pool(gb_per_process=gb_per_process, num_processes=num_processes)

    for ccd in xmm.ccds:
        name = "ccd{:02d}".format(ccd.num)
        xml_file_path = ccd.file_location
        run_sixte(name, simput_file_path, xml_file_path, rundir, exposure, ra, dec, rollangle, verbose)

        # pool.apply_async(run_sixte, args=(name, simput_file_path, xml_file_path, rundir, exposure, ra, dec, rollangle,))

    # # Close the pool
    # pool.close()
    # pool.join()
