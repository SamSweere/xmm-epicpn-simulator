import os
import time
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List, Union, Tuple
from uuid import uuid4
from warnings import warn

import numpy as np

from src.simput.agn import get_fluxes
from src.simput.gen import background, generate_exposure_map, simput_ps
from src.simput.gen.image import simput_image
from src.simput.utils import get_spectrumfile
from src.simput.utils import merge_simputs
from utils.file_utils import compress_gzip
from utils.log import elog


_spectrum_ds_file = Path("/home/bojantodorkov/Projects/xmm-epicpn-simulator/res/simput/pntffg_spectrum.ds")


def create_background(
        run_dir: Path,
        num: int = 1,
        verbose: bool = True
) -> List[Path]:
    output_files = []

    # TODO Generation does not work
    for _ in range(num):
        output_file = background(run_dir=run_dir,
                                 spectrum_file=_spectrum_ds_file,
                                 emin=0.15,
                                 emax=15.0,
                                 verbose=verbose)
        output_files.append(output_file)

    return output_files


def create_exposure_map(
        run_dir,
        num: int = 1,
        verbose: bool = True
) -> List[Path]:
    output_files = []

    for _ in range(num):
        output_file = generate_exposure_map(run_dir=run_dir,
                                            spectrum_file=_spectrum_ds_file,
                                            emin=0.15,
                                            emax=15.0,
                                            verbose=verbose)
        output_files.append(output_file)

    return output_files


def create_random_sources(
        run_dir: Path,
        num: int = 1,
        num_sources: int = 10,
        keep_files: bool = False,
        verbose: bool = True
) -> List[Path]:
    """
    Generates an image with multiple random sources in it.
    """
    output_files = []
    for _ in range(num):
        simput_files = []

        spectrum_file = get_spectrumfile(run_dir=run_dir, norm=0.01)

        # emin = 0.15, emax = 15.0 is ideal
        # For some reason if I go into the emin=0.5 and emax=2.0 sixte will not complete anymore
        emin = 0.15
        emax = 15.0

        name = f"random_n_{num_sources}_{emin}eV_p1_{emax}eV"
        output_file = run_dir / f"{name}.simput"

        # Create the point sources
        for i in range(num_sources):
            print(f"Generating point source {i}/{num_sources - 1}")
            simput_file = run_dir / f"ps_{i}.simput"
            simput_ps(emin=emin,
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
                                      keep_files=keep_files,
                                      verbose=verbose)

        output_files.append(random_source)

    return output_files


def create_agn_sources(
        run_dir: Path,
        num: int = 1,
        keep_files: bool = False,
        verbose: bool = True
):
    if not isinstance(num, int):
        raise ValueError
    output_files = []
    emin = 0.5
    emax = 2.0

    # Get the fluxes from the agn distribution
    fluxes = get_fluxes()

    for _ in range(num):
        # Use the current time as id, such that clashes don't happen
        unique_id = uuid4().int
        output_file_path = run_dir / f"agn_{unique_id}_p0_{emin}ev_p1_{emax}ev.simput"
        simput_files = []
        # Get the spectrum file
        spectrum_file = get_spectrumfile(run_dir=run_dir, norm=0.001, verbose=verbose)
        for i, flux in enumerate(fluxes):
            if verbose:
                print("Creating source with flux:", flux)
            output_file = run_dir / f"ps_{unique_id}_{i}.simput"
            output_file = simput_ps(emin=emin,
                                    emax=emax,
                                    output_file=output_file,
                                    src_flux=flux,
                                    xspec_file=spectrum_file,
                                    offset='random',
                                    verbose=verbose)
            simput_files.append(output_file)
        output_file = merge_simputs(simput_files=simput_files,
                                    output_file=output_file_path,
                                    keep_files=keep_files,
                                    verbose=verbose)
        output_files.append(output_file)

    return output_files


def create_test_grid(
        run_dir: Path,
        num: int = 1,
        flux: float = 1.0e-13,
        step_size: int = 10,
        center_offset: Tuple[float, float] = (0.005, 0.015),
        keep_files: bool = False,
        verbose=True
) -> List[Path]:
    output_files = []
    # emin = 0.15, emax = 15.0 is ideal
    # For some reason if I go into the emin=0.5 and emax=2.0 sixte will not complete anymore
    emin = 0.15
    emax = 15.0

    xmm_fov = 0.5
    x_loc = np.linspace(-xmm_fov / 2, xmm_fov / 2, step_size)
    y_loc = np.linspace(-xmm_fov / 2, xmm_fov / 2, step_size)
    for _ in range(num):
        unique_id = uuid4().int
        name = f"test_grid_f_{flux}_step_{step_size}_{unique_id}.simput"
        output_file = run_dir / name

        counter = 0
        simput_files = []
        spectrum_file = get_spectrumfile(run_dir=run_dir, norm=0.01, verbose=verbose)
        # Create one point source at (0, 0)
        if verbose:
            print(f"Generating point source at (0, 0)")
        simput_file_path = run_dir / f"ps_center_{unique_id}.simput"
        simput_ps(emin=emin, emax=emax, center_point=(0.0, 0.0), xspec_file=spectrum_file,
                  output_file=simput_file_path, src_flux=flux, offset=(0.0, 0.0), verbose=verbose)
        simput_files.append(simput_file_path)

        # Create the point sources
        for x in x_loc:
            for y in y_loc:
                if verbose:
                    print(f"Generating point source {counter}/{step_size ** 2}")
                ps_simput_file_name = f"ps_{counter}_{unique_id}.simput"
                simput_file_path = run_dir / ps_simput_file_name
                simput_ps(emin=emin, emax=emax, center_point=center_offset, xspec_file=spectrum_file,
                          output_file=simput_file_path, src_flux=flux, offset=(x, y), verbose=verbose)
                simput_files.append(simput_file_path)
                counter += 1

        output_file = merge_simputs(simput_files=simput_files,
                                    output_file=output_file,
                                    keep_files=keep_files,
                                    verbose=verbose)

        output_files.append(output_file)

    return output_files


def simput_generate(
        mode: str,
        img_settings: Union[dict, int],
        tmp_dir: Path,
        output_dir: Path,
        keep_files: bool = False,
        verbose: bool = True
) -> None:
    if mode != "test_grid":
        return
    with TemporaryDirectory(dir=tmp_dir) as temp:
        run_dir = Path(temp)

        if mode == "test_grid":
            file_names = create_test_grid(run_dir=run_dir,
                                          num=img_settings,
                                          keep_files=keep_files,
                                          verbose=verbose)
        elif mode == "agn":
            file_names = create_agn_sources(run_dir=run_dir,
                                            num=img_settings,
                                            keep_files=keep_files,
                                            verbose=verbose)
        elif mode == "img":
            file_names = simput_image(run_dir=run_dir,
                                      img_settings=img_settings,
                                      keep_files=keep_files,
                                      verbose=verbose)
        elif mode == "background":
            # TODO
            file_names = create_background(run_dir=run_dir,
                                           num=img_settings,
                                           verbose=verbose)
            pass
        elif mode == "exposure_map":
            # TODO
            # file_names = create_exposure_map(run_dir=run_dir,
            #                                  num=img_settings,
            #                                  verbose=verbose)
            pass
        elif mode == "random":
            file_names = create_random_sources(run_dir=run_dir,
                                               num=img_settings,
                                               keep_files=keep_files,
                                               verbose=verbose)
        else:
            e_message = f"Mode {mode} is not supported"
            elog(e_message)
            raise ValueError(e_message)

        for file_name in file_names:
            # Compress the simput file and move it to the correct output dir
            compressed_file = output_dir / f"{file_name.name}.gz"
            if compressed_file.exists():
                m = f"WARNING: Simput file {compressed_file.resolve()} already exists, skipping"
                warn(m)
                # plog(m, verbose)
            else:
                compress_gzip(in_file_path=file_name, out_file_path=compressed_file)

            # plog(f"Saved compressed simput file at: {compressed_file}", verbose)
