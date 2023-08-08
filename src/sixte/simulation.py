from pathlib import Path

from src.utils.external_run import run_headas_command


def run_sixte(
        name: str,
        simput_file_path,
        xml_file_path,
        output_dir: Path,
        exposure,
        ra,
        dec,
        rollangle,
        verbose: bool
):
    if verbose:
        print("--------------------")
        print("Running " + name)

    raw_filepath = output_dir / f"{name}_raw.fits"
    evt_filepath = output_dir / f"{name}_evt.fits"

    sixte_command = f"runsixt RawData={raw_filepath} EvtFile={evt_filepath} " \
                    f"XMLFile={xml_file_path} Simput={simput_file_path} Exposure={exposure} RA={ra} " \
                    f"Dec={dec} rollangle={rollangle} clobber=yes"

    run_headas_command(sixte_command, verbose=verbose)


def run_sixte_all_ccds(
        xmm,
        simput_file_path,
        rundir,
        exposure,
        ra=0.0,
        dec=0.0,
        rollangle=0.0,
        verbose=True
):
    # Run the sixte in parallel
    # TODO Enable multiprocessing
    # pool = get_multiprocessing_pool(gb_per_process=gb_per_process, num_processes=num_processes)

    for ccd in xmm.ccds:
        name = "ccd{:02d}".format(ccd.num)
        # xml_file_path = ccd.file_location
        xml_file_path = "/home/bojantodorkov/simput/share/sixte/instruments/xmm/epicpn/fullframe_thinfilter_ccd00.xml"
        run_sixte(name, simput_file_path, xml_file_path, rundir, exposure, ra, dec, rollangle, verbose)

        # pool.apply_async(run_sixte, args=(name, simput_file_path, xml_file_path, rundir, exposure, ra, dec, rollangle,))

    # # Close the pool
    # pool.close()
    # pool.join()
