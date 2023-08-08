# Simulate two point sources with power-law spectra
from pathlib import Path
from typing import Tuple, Union

import numpy as np

from src.utils.external_run import run_headas_command


def simput_ps(
        emin: float,
        emax: float,
        output_file: Path,
        xspec_file: Path,
        center_point: Tuple[float, float] = (0.0, 0.0),
        src_flux: Union[float, str] = 1.0e-12,
        offset: Union[Tuple[float, float], str] = (0.0, 0.0),
        verbose: bool = True
) -> Path:
    """
    Generates a single point-source
    """
    if isinstance(offset, str) and offset != "random":
        raise ValueError(f"Value of offset is unknown string \"{offset}\"!")

    if isinstance(src_flux, str) and src_flux != "random":
        raise ValueError(f"Value of src_flux is unknown string \"{src_flux}\"!")

    rng = np.random.default_rng()
    # TODO these values are only based on the tutorial values, no thought if they are realistic
    if src_flux == 'random':
        src_flux = rng.uniform(low=1.0e-13, high=1.0e-10)

    # Epic pn has a fov of 30 arcmin = 0.5 degrees.
    pn_fov = 0.5

    # Randomly position the point source within the fov
    if offset == 'random':
        offset = rng.uniform(low=-1.0 * pn_fov / 2, high=pn_fov / 2, size=2)

    location = (center_point[0] + offset[0], center_point[1] + offset[1])
    ra = location[0]
    if ra < 0:
        ra = 360 + ra
    dec = location[1]
    # if dec < 0:
    #     dec = 90 + dec

    simput_command = f"simputfile RA={ra} Dec={dec} XSPECFile={xspec_file.resolve()} Emin={emin} Emax={emax} " \
                     f"srcFlux={src_flux} Simput={output_file.resolve()}"

    run_headas_command(simput_command, verbose=verbose)
    return output_file
