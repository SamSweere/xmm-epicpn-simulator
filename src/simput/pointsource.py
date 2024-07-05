from pathlib import Path

import numpy as np

from src.sixte import commands
from src.xmm.tools import get_fov


def create_pointsource(
    emin: float,
    emax: float,
    output_file: Path,
    xspec_file: Path,
    center_point: tuple[float, float] = (0.0, 0.0),
    src_flux: float | str = 1.0e-12,
    offset: tuple[float, float] | str = (0.0, 0.0),
) -> Path:
    """
    Generates a single point-source

    Returns:
        Path: Path to the file containing the generated single point-source
    """
    if isinstance(offset, str) and offset != "random":
        raise ValueError(f'Value of offset is unknown string "{offset}"!')

    if isinstance(src_flux, str) and src_flux != "random":
        raise ValueError(f'Value of src_flux is unknown string "{src_flux}"!')

    rng = np.random.default_rng()
    # TODO these values are only based on the tutorial values, no thought if they are realistic
    if isinstance(src_flux, str) and src_flux == "random":
        src_flux = rng.uniform(low=1.0e-13, high=1.0e-10)

    # Randomly position the point source within the fov
    if isinstance(offset, str) and offset == "random":
        # The FOV is the same for EPN, EMOS1, and EMOS2
        fov = get_fov("epn")
        offset = rng.uniform(low=-1.0 * fov / 2, high=fov / 2, size=2)

    location = (center_point[0] + offset[0], center_point[1] + offset[1])
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
