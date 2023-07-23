# Simulate two point sources with power-law spectra
import numpy as np

from src.utils.external_run import run_headas_command


def generate_point_source(emin, emax, simput_file_path, xspec_file, center_point=(0.0, 0.0), src_flux=1.0e-12,
                          offset=(0.0, 0.0), verbose=True):
    rng = np.random.default_rng()
    # If the src_flux of offset is None then the values are randomly generated
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

    # We need the xspec from headas
    simput_command = f"simputfile RA={ra} Dec={dec} XSPECFile={xspec_file} Emin={emin} Emax={emax} " \
                     f"srcFlux={src_flux} Simput={simput_file_path} "

    run_headas_command(simput_command, verbose=verbose)
