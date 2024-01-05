import os
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

from yt.utilities.logger import set_log_level
from yt.loaders import load
from yt.fields.xray_emission_fields import add_xray_emissivity_field
from yt.visualization.fits_image import FITSOffAxisProjection, FITSSlice
from loguru import logger

set_log_level(
    50
)  # Turn of logging to console since it could lead to a deadlock when using multiprocessing


def cutout_to_xray_fits(
    cutout: Path,
    output_dir: Path,
    sub,
    mode_dict: Dict[str, list],
    cloudy_emissivity_root: Path,
    widths: List[Tuple[float, str]],
    resolutions: List[int],
    emin: float,
    emax: float,
    redshift: float,
    overwrite: bool = True,
    fail_on_error: bool = False,
    consume_data: bool = False,
):
    fits_filename_prefix = cutout.stem
    output_dir = (
        output_dir / cutout.parts[-3] / cutout.parts[-2]
    )  # Save in output_dir / TNG set / snapshot num
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.debug(
        f"Processing cutout: {cutout.resolve()}. Arguments: mode={mode_dict}; emin={emin}; emax={emax}; "
        f"widths={widths}; resolution={resolutions}; redshift={redshift}; overwrite={overwrite}; "
        f"output_dir={output_dir}; cutout_datafolder={cutout.parent}"
    )
    start = datetime.now()

    # Force turn off all logging from yt via redirecting stdout and stderr into /dev/null
    with open(os.devnull, "w") as f, redirect_stdout(f), redirect_stderr(f):
        ds = load(f"{cutout.resolve()}", default_species_fields="ionized")

        add_xray_emissivity_field(
            ds, emin, emax, redshift, data_dir=cloudy_emissivity_root
        )

        for mode in mode_dict.keys():
            axes = mode_dict[mode]
            if not axes:
                continue

            for axis in axes:
                for w in widths:
                    for r in resolutions:
                        logger.debug(f"mode={mode}, axis={axis}, w={w}, r={r}")
                        fits_filename_suffix = (
                            f"_m_{mode}_r_{r}_w_{w[0]}{w[1]}_n_{axis}"
                        )
                        fits_filename = (
                            f"{fits_filename_prefix}{fits_filename_suffix}.fits"
                        )
                        fits_path = output_dir / fits_filename
                        fits_path = fits_path.resolve()

                        if fits_path.exists() and not overwrite:
                            raise FileExistsError(
                                f"{fits_path} already exists and `overwrite` is False"
                            )
                        try:
                            if mode == "proj":
                                yt_fits = FITSOffAxisProjection(
                                    ds,
                                    normal=w,
                                    fields=(
                                        "gas",
                                        f"xray_photon_intensity_{emin}_{emax}_keV",
                                    ),
                                    center=(sub["pos_x"], sub["pos_y"], sub["pos_z"]),
                                    width=w,
                                    image_res=r,
                                )
                            else:
                                yt_fits = FITSSlice(
                                    ds,
                                    axis=axis,
                                    fields=(
                                        "gas",
                                        f"xray_photon_intensity_{emin}_{emax}_keV",
                                    ),
                                    center=(sub["pos_x"], sub["pos_y"], sub["pos_z"]),
                                    width=w,
                                    image_res=r,
                                )

                            yt_fits.update_header(field="all", key="AXIS", value=axis)
                            yt_fits.update_header(field="all", key="WIDTH", value=w)
                            yt_fits.update_header(
                                field="all", key="REDSHIFT", value=redshift
                            )
                            yt_fits.update_header(field="all", key="EMIN", value=emin)
                            yt_fits.update_header(field="all", key="EMAX", value=emax)
                            yt_fits.writeto(f"{fits_path}", overwrite=overwrite)
                        except Exception as e:
                            if fail_on_error:
                                logger.exception(
                                    f"Failed to process {cutout.name} with error:\n"
                                    f"{e}"
                                )
                                raise
                            else:
                                logger.warning(
                                    f"Failed to process {cutout.name} with error:\n"
                                    f"{e}"
                                )
    if consume_data:
        cutout.unlink()
    end = datetime.now()
    logger.info(f"Processing of {cutout} took {end - start}.")
