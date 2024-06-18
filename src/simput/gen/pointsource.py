from pathlib import Path

import numpy as np

from src.sixte import commands
from src.xmm.utils import get_fov

import warnings


def simput_ps(
    emin: float,
    emax: float,
    output_file: Path,
    xspec_file: Path,
    # center_point: tuple[float, float] = (0.0, 0.0),
    location: tuple[float, float],
    src_flux: float | str = 1.0e-12,
    verbose: bool = True,
) -> Path:
    """
    Generates a single point-source

    Returns:
        Path: Path to the file containing the generated single point-source
    """
   
    if isinstance(src_flux, str) and src_flux != "random":
        raise ValueError(f'Value of src_flux is unknown string "{src_flux}"!')

    rng = np.random.default_rng()
    # TODO these values are only based on the tutorial values, no thought if they are realistic
    # TODO: make sure that this also works for deblending, right now it only works if non-random
    if src_flux == "random":
        src_flux = rng.uniform(low=1.0e-13, high=1.0e-10)
        warnings.warn("The deblending option is not properly implemented for random fluxes!", UserWarning)
        

    ra = location[0]
    if ra < 0:
        ra = 360 + ra

    dec = location[1]
    # if dec < 0:
    #     dec = 90 + dec

    commands.simputfile(
        simput=output_file, ra=ra, dec=dec, src_flux=src_flux, emin=emin, emax=emax, xspec_file=xspec_file
    )

    return output_file.resolve()
