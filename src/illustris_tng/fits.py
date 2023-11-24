import os
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import yt
from loguru import logger

yt.set_log_level(50)  # Turn of logging to console since it could lead to a deadlock when using multiprocessing


def cutout_to_xray_fits(
        cutout: Path,
        output_dir: Path,
        sub,
        mode_dict: Dict[str, list],
        cloudy_emissivity_root: Path,
        emin=0.5,
        emax=2.0,
        width: List[list] = None,
        resolution: List[int] = None,
        redshift=0.05,
        overwrite=False
):
    if width is None:
        width = [(1.0, "unitary")]
    if resolution is None:
        resolution = [2048]
    fits_filename_prefix = cutout.stem
    output_dir = output_dir / cutout.parts[-3] / cutout.parts[-2]  # Save in output_dir / TNG set / snapshot num
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Processing cutout: {cutout.resolve()}. Arguments: mode={mode_dict}; emin={emin}; emax={emax}; "
                f"width={width}; resolution={resolution}; redshift={redshift}; overwrite={overwrite}; "
                f"output_dir={output_dir}; cutout_datafolder={cutout.parent}")
    start = datetime.now()

    # Force turn off all logging from yt via redirecting stdout and stderr into /dev/null
    with open(os.devnull, "w") as f, redirect_stdout(f), redirect_stderr(f):
        ds = yt.load(f"{cutout.resolve()}", default_species_fields="ionized")

        yt.add_xray_emissivity_field(ds, emin, emax, redshift, data_dir=cloudy_emissivity_root)

        for mode, axes in mode_dict.items():
            for axis in axes:
                for w in width:
                    w = tuple(w)
                    for r in resolution:
                        fits_filename_suffix = f"_m_{mode}_r_{r}_w_{w[0]}{w[1]}_n_{axis}"
                        fits_filename = f"{fits_filename_prefix}{fits_filename_suffix}.fits"
                        fits_path = output_dir / fits_filename
                        fits_path = fits_path.resolve()

                        if fits_path.exists() and not overwrite:
                            raise FileExistsError(f"{fits_path} already exists and `overwrite` is False")
                        try:
                            if mode == "proj":
                                normal = []
                                if axis == "x":
                                    normal = [1, 0, 0]
                                if axis == "y":
                                    normal = [0, 1, 0]
                                if axis == "z":
                                    normal = [0, 0, 1]
                                yt_fits = yt.FITSOffAxisProjection(ds, normal=normal,
                                                                   fields=('gas',
                                                                           f'xray_photon_intensity_{emin}_{emax}_keV'),
                                                                   center=(sub['pos_x'],
                                                                           sub['pos_y'],
                                                                           sub['pos_z']),
                                                                   width=w, image_res=r)
                            else:
                                yt_fits = yt.FITSSlice(ds, axis=axis,
                                                       fields=('gas', f'xray_photon_intensity_{emin}_{emax}_keV'),
                                                       center=(sub['pos_x'], sub['pos_y'], sub['pos_z']),
                                                       width=w, image_res=r)
                        except Exception as e:
                            logger.exception(f"ERROR\tFailed to process {cutout.name} with error:\n"
                                             f"{e}")

                        yt_fits.update_header(field="all", key="AXIS", value=axis)
                        yt_fits.update_header(field="all", key="WIDTH", value=w)
                        yt_fits.update_header(field="all", key="REDSHIFT", value=redshift)
                        yt_fits.update_header(field="all", key="EMIN", value=emin)
                        yt_fits.update_header(field="all", key="EMAX", value=emax)
                        yt_fits.writeto(f"{fits_path}", overwrite=overwrite)
    end = datetime.now()
    logger.info(f"DONE\tProcessing of {cutout} took {end - start}.")
