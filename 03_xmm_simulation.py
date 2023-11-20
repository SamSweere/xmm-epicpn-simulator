import json
from argparse import ArgumentParser
from functools import partial
from multiprocessing.pool import Pool
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict, Union

from loguru import logger

from src.simput.utils import get_simputs
from src.sixte.simulator import run_xmm_simulation
from src.xmm.utils import create_psf_file, create_vinget_file, create_xml_files
from src.xmm_utils.multiprocessing import get_num_processes
from src.xmm_utils.run_utils import configure_logger, handle_error
from datetime import timedelta

logger.remove()


def run(path_to_cfg: Union[Path, Dict[str, dict]]) -> None:
    if isinstance(path_to_cfg, Path):
        with open(path_to_cfg, "r") as f:
            cfg: Dict[str, dict] = json.load(f)

    env_cfg = cfg["environment"]
    mp_cfg = cfg["multiprocessing"]
    instrument_cfg = cfg["instrument"]

    log_dir = Path(env_cfg["log_directory"]).expanduser()

    debug = env_cfg["debug"]

    configure_logger(log_dir=log_dir, log_name="03_xmm_simulation.log", enqueue=True, debug=debug,
                     rotation=timedelta(hours=1))

    if not debug:
        logger.info(f"Since 'debug' is set to 'false' the simulation will be run asynchronously.")

    working_directory = Path(env_cfg["working_directory"]).expanduser()
    working_directory.mkdir(exist_ok=True, parents=True)
    out_dir = working_directory / "xmm_sim_dataset"
    out_dir.mkdir(exist_ok=True)
    simput_path = working_directory / "simput"
    if not simput_path.exists():
        raise NotADirectoryError(f"Simput directory '{simput_path.resolve()}' not found! "
                                 f"Please make sure that the simput files have been generated and that "
                                 f"the working directory is set correctly.")

    instrument_names = instrument_cfg["instrument_names"]
    res_mults = instrument_cfg["res_mults"]

    with TemporaryDirectory(prefix="xml") as xml_dir, TemporaryDirectory(prefix="sim") as sim_dir:
        xml_dir = Path(xml_dir)
        sim_dir = Path(sim_dir)

        # Create all needed directories
        for instrument_name in instrument_names:
            instrument_dir = xml_dir / instrument_name
            for res_mult in res_mults:
                res_mult_dir = instrument_dir / instrument_cfg["filter"] / f"{res_mult}x"
                res_mult_dir.mkdir(exist_ok=True, parents=True)

        with Pool(len(instrument_names)) as pool:
            mp_apply = pool.apply if debug else partial(pool.apply_async, error_callback=handle_error)
            logger.info("START\tCreating all PSF files.")
            for instrument_name in instrument_names:
                for res_mult in res_mults:
                    arguments = (xml_dir, instrument_name, res_mult)
                    mp_apply(create_psf_file, arguments)
            pool.close()
            pool.join()
            logger.info("DONE\tPSF files have been created.")

        with Pool(len(instrument_names)) as pool:
            mp_apply = pool.apply if debug else partial(pool.apply_async, error_callback=handle_error)
            logger.info("START\tCreating all vignetting files.")
            for instrument_name in instrument_names:
                arguments = (xml_dir, instrument_name)
                mp_apply(create_vinget_file, arguments)
            pool.close()
            pool.join()
            logger.info("DONE\tVignetting files have been created.")

        logger.info("START\tCreating all XML files.")
        for instrument_name in instrument_names:
            for res_mult in res_mults:
                create_xml_files(xml_dir, instrument_name, res_mult, instrument_cfg["filter"],
                                 instrument_cfg["sim_separate_ccds"], instrument_cfg["wait_time"])
        logger.info("DONE\tXML files have been created.")

        mode_amount_dict: Dict[str, int] = instrument_cfg["mode"]

        max_exposure = instrument_cfg['max_exposure']
        max_exp_str = f"{round(max_exposure / 1000)}ks"
        xmm_filter = instrument_cfg["filter"]

        with Pool(get_num_processes(mp_conf=mp_cfg)) as pool:
            mp_apply = pool.apply if debug else partial(pool.apply_async, error_callback=handle_error)
            for instrument_name in instrument_names:
                mode_simput_files = get_simputs(instrument_name=instrument_name, simput_path=simput_path,
                                                mode_amount_dict=mode_amount_dict, order=instrument_cfg["order"])
                if "background" in mode_simput_files:
                    # Since we only have one background but want multiple simulations of it repeat it
                    mode_simput_files["background"] = mode_simput_files["background"] * mode_amount_dict["background"]

                xmm_filter_dir = out_dir / instrument_name / xmm_filter
                max_exp_dir = xmm_filter_dir / max_exp_str
                for res_mult in res_mults:
                    res_str = f"{res_mult}x"
                    for mode, files in mode_simput_files.items():
                        res_dir = max_exp_dir / mode / res_str
                        for file in files:
                            img_name = f"{file.name.replace('.simput.gz', '')}_mult_{res_mult}"
                            final_compressed_img_name = f"{img_name}_{max_exp_str}_p_0-0.fits.gz"
                            final_compressed_file_path = res_dir / final_compressed_img_name
                            if not debug:
                                if mode != "background":
                                    # Check if the max exposure time already exists,
                                    # if so we know that this one was already generated,
                                    # and we can skip its generation
                                    if final_compressed_file_path.exists():
                                        continue
                            simulation_args = (xml_dir.resolve(), file.resolve(), img_name, mode, sim_dir.resolve(),
                                               xmm_filter_dir.resolve(), res_mult, max_exposure, instrument_name,
                                               xmm_filter, instrument_cfg["sim_separate_ccds"], debug,
                                               env_cfg["verbose"])

                            mp_apply(run_xmm_simulation, simulation_args)
            pool.close()
            pool.join()
            logger.info("Done!")


if __name__ == '__main__':
    parser = ArgumentParser(prog="", description="")
    parser.add_argument("-p", "--config_path", type=Path, required=True, help="Path to config file.")

    args = parser.parse_args()
    run(path_to_cfg=args.config_path)
