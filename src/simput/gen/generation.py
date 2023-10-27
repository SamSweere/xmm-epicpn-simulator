import itertools
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List, Literal
from uuid import uuid4

import numpy as np
from loguru import logger

from src.simput import constants
from src.simput.agn import get_fluxes
from src.simput.gen import background, exposure_map, simput_ps
from src.simput.gen.image import simput_image
from src.simput.utils import get_spectrumfile, merge_simputs
from src.xmm.utils import get_fov_for_instrument
from src.xmm_utils.file_utils import compress_gzip


def create_background(
        instrument_name: Literal["epn", "emos1", "emos2"],
        emin: float,
        emax: float,
        run_dir: Path,
        spectrum_file: Path,
        num: int = 1,
        verbose: bool = True
) -> List[Path]:
    if num == 1:
        output_files = [background(run_dir=run_dir,
                                   spectrum_file=spectrum_file,
                                   instrument_name=instrument_name,
                                   emin=emin,
                                   emax=emax,
                                   verbose=verbose)]
    else:
        output_files = [background(run_dir=run_dir,
                                   spectrum_file=spectrum_file,
                                   instrument_name=instrument_name,
                                   emin=emin,
                                   emax=emax,
                                   suffix=i,
                                   verbose=verbose) for i in range(num)]

    return output_files


def create_exposure_map(
        instrument_name: Literal["epn", "emos1", "emos2"],
        emin: float,
        emax: float,
        run_dir: Path,
        spectrum_file: Path,
        num: int = 1,
        verbose: bool = True
) -> List[Path]:
    if num == 1:
        output_files = [exposure_map(run_dir=run_dir,
                                     spectrum_file=spectrum_file,
                                     instrument_name=instrument_name,
                                     emin=emin,
                                     emax=emax,
                                     verbose=verbose)]
    else:
        output_files = [exposure_map(run_dir=run_dir,
                                     spectrum_file=spectrum_file,
                                     instrument_name=instrument_name,
                                     emin=emin,
                                     emax=emax,
                                     suffix=i,
                                     verbose=verbose) for i in range(num)]

    return output_files


def create_random_sources(
        instrument_name: Literal["epn", "emos1", "emos2"],
        emin: float,
        emax: float,
        run_dir: Path,
        num: int = 1,
        num_sources: int = 10,
        verbose: bool = True
) -> List[Path]:
    """
    Generates an image with multiple random sources in it.
    """
    spectrum_file = get_spectrumfile(run_dir=run_dir, verbose=verbose)

    output_files = []
    for _ in range(num):
        simput_files = []

        name = f"random_n_{num_sources}_{emin}eV_p1_{emax}eV"
        output_file = run_dir / f"{name}.simput"

        # Create the point sources
        for i in range(num_sources):
            simput_file = run_dir / f"ps_{i}.simput"
            simput_ps(instrument_name=instrument_name,
                      emin=emin,
                      emax=emax,
                      center_point=(0, 0),
                      xspec_file=spectrum_file,
                      output_file=simput_file,
                      src_flux='random',
                      offset='random',
                      verbose=verbose)
            simput_files.append(simput_file)

        random_source = merge_simputs(simput_files=simput_files,
                                      output_file=output_file,
                                      verbose=verbose)

        output_files.append(random_source)

    return output_files


def create_agn_sources(
        instrument_name: Literal["epn", "emos1", "emos2"],
        emin: float,
        emax: float,
        run_dir: Path,
        agn_counts_file: Path,
        num: int = 1,
        verbose: bool = True
):
    output_files = []

    # Get the spectrum file
    spectrum_file = get_spectrumfile(run_dir=run_dir, norm=0.001, verbose=verbose)
    # Get the fluxes from the agn distribution
    fluxes = get_fluxes(agn_counts_file)

    for _ in range(num):
        # Use the current time as id, such that clashes don't happen
        unique_id = uuid4().int
        output_file_path = run_dir / f"agn_{unique_id}_p0_{emin}ev_p1_{emax}ev.simput"
        simput_files = []

        for i, flux in enumerate(fluxes):
            if verbose:
                logger.info(f"Creating AGN with flux={flux}")
            output_file = run_dir / f"ps_{unique_id}_{i}.simput"
            output_file = simput_ps(instrument_name=instrument_name,
                                    emin=emin,
                                    emax=emax,
                                    output_file=output_file,
                                    src_flux=flux,
                                    xspec_file=spectrum_file,
                                    offset='random',
                                    verbose=verbose)
            simput_files.append(output_file)
        output_file = merge_simputs(simput_files=simput_files,
                                    output_file=output_file_path,
                                    verbose=verbose)
        output_files.append(output_file)

    return output_files


def create_test_grid(
        instrument_name: Literal["epn", "emos1", "emos2"],
        emin: float,
        emax: float,
        run_dir: Path,
        num: int = 1,
        flux: float = 1.0e-13,
        step_size: int = 10,
        verbose=True
) -> List[Path]:
    spectrum_file = get_spectrumfile(run_dir=run_dir, verbose=verbose)

    output_files = []

    fov = get_fov_for_instrument(instrument_name=instrument_name)
    loc = list(itertools.product(np.linspace(-fov / 2, fov / 2, step_size), repeat=2))
    loc.append((0.0, 0.0))
    for _ in range(num):
        unique_id = uuid4().int
        name = f"test_grid_f_{flux}_step_{step_size}_{unique_id}.simput"
        output_file = run_dir / name

        simput_files = []
        for x, y in loc:
            unique_id = uuid4().int
            if verbose:
                logger.info(f"Generating point source at x={x}, y={y}.")
            ps_simput_file_name = f"ps_{x}_{y}_{unique_id}.simput"
            simput_file_path = run_dir / ps_simput_file_name
            simput_ps(instrument_name=instrument_name,
                      emin=emin,
                      emax=emax,
                      xspec_file=spectrum_file,
                      output_file=simput_file_path,
                      src_flux=flux,
                      offset=(x, y),
                      verbose=verbose)
            simput_files.append(simput_file_path)

        output_file = merge_simputs(simput_files=simput_files,
                                    output_file=output_file,
                                    verbose=verbose)

        output_files.append(output_file)

    return output_files


def simput_generate(
        instrument_name: Literal["epn", "emos1", "emos2"],
        emin: float,
        emax: float,
        mode: str,
        img_settings: dict,
        tmp_dir: Path,
        output_dir: Path,
        verbose: bool = True
) -> None:
    if mode not in constants.available_modes:
        raise ValueError(f"Unkown mode '{mode}'! Available modes: {constants.available_modes}.")

    with TemporaryDirectory(dir=tmp_dir) as temp:
        run_dir = Path(temp)

        if mode == "test_grid":
            file_names = create_test_grid(instrument_name=instrument_name,
                                          emin=emin,
                                          emax=emax,
                                          run_dir=run_dir,
                                          num=img_settings["num"],
                                          verbose=verbose)

        if mode == "agn":
            if img_settings["agn_counts_file"] is None:
                raise FileNotFoundError(f"Path to agn_counts.cgi cannot be None!")
            file_names = create_agn_sources(instrument_name=instrument_name,
                                            emin=emin,
                                            emax=emax,
                                            run_dir=run_dir,
                                            agn_counts_file=img_settings["agn_counts_file"],
                                            num=img_settings["num"],
                                            verbose=verbose)

        if mode == "img":
            file_names = simput_image(instrument_name=instrument_name,
                                      emin=emin,
                                      emax=emax,
                                      run_dir=run_dir,
                                      img_settings=img_settings,
                                      verbose=verbose)

        if mode == "background":
            if img_settings["spectrum_file"] is None:
                raise FileNotFoundError(f"Path to pnttfg_spectrum.ds cannot be None!")
            file_names = create_background(instrument_name=instrument_name,
                                           emin=emin,
                                           emax=emax,
                                           run_dir=run_dir,
                                           spectrum_file=img_settings["spectrum_file"],
                                           num=img_settings["num"],
                                           verbose=verbose)

        if mode == "exposure_map":
            if img_settings["spectrum_file"] is None:
                raise FileNotFoundError(f"Path to pnttfg_spectrum.ds cannot be None!")
            file_names = create_exposure_map(instrument_name=instrument_name,
                                             emin=emin,
                                             emax=emax,
                                             run_dir=run_dir,
                                             spectrum_file=img_settings["spectrum_file"],
                                             num=img_settings["num"],
                                             verbose=verbose)

        if mode == "random":
            file_names = create_random_sources(instrument_name=instrument_name,
                                               emin=emin,
                                               emax=emax,
                                               run_dir=run_dir,
                                               num=img_settings["num"],
                                               verbose=verbose)

        for file_name in file_names:
            # Compress the simput file and move it to the correct output dir
            compressed_file = output_dir / f"{file_name.name}.gz"
            if compressed_file.exists():
                logger.warning(f"Simput file {compressed_file.resolve()} already exists, skipping.")
            else:
                compress_gzip(in_file_path=file_name, out_file_path=compressed_file)
