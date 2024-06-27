import os
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime
from pathlib import Path

from loguru import logger
from yt.fields.xray_emission_fields import add_xray_emissivity_field
from yt.loaders import load
from yt.utilities.logger import set_log_level
from yt.visualization.fits_image import FITSOffAxisProjection, FITSSlice

from src.config import EnergyCfg, EnvironmentCfg
from src.tools.files import compress_gzip

set_log_level(50)


def cutout_to_xray_fits(
    cutout: Path,
    output_dir: Path,
    sub,
    mode_dict: dict[str, list],
    cloudy_emissivity_root: Path,
    widths: list[tuple[float, str]],
    resolutions: list[int],
    energies: EnergyCfg,
    redshift: float,
    environment: EnvironmentCfg,
) -> list[Path]:
    fits_filename_prefix = cutout.stem
    output_dir = output_dir / cutout.parts[-3] / cutout.parts[-2]  # Save in output_dir / TNG set / snapshot num
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.debug(
        f"Processing cutout: {cutout.resolve()}. Arguments: mode={mode_dict}; emin={energies.emin}; "
        f"emax={energies.emax}; widths={widths}; resolution={resolutions}; redshift={redshift}; "
        f"overwrite={environment.overwrite}; output_dir={output_dir}; cutout_datafolder={cutout.parent}"
    )
    fits = []
    start = datetime.now()

    # Force turn off all logging from yt via redirecting stdout and stderr into /dev/null
    with open(os.devnull, "w") as f, redirect_stdout(f), redirect_stderr(f):
        ds = load(cutout, default_species_fields="ionized")

        add_xray_emissivity_field(ds, energies.emin, energies.emax, redshift, data_dir=cloudy_emissivity_root)

        for mode in mode_dict:
            axes = mode_dict[mode]
            if not axes:
                continue

            for axis in axes:
                for w in widths:
                    for r in resolutions:
                        logger.debug(f"mode={mode}, axis={axis}, w={w}, r={r}")
                        axis_str = axis if isinstance(axis, str) else f"{axis[0]}_{axis[1]}_{axis[2]}"
                        fits_filename_suffix = f"_m_{mode}_r_{r}_w_{w[0]}{w[1]}_n_{axis_str}"
                        fits_filename = f"{fits_filename_prefix}{fits_filename_suffix}.fits"
                        fits_path = output_dir / fits_filename
                        compressed_path = fits_path.with_suffix(".fits.gz")

                        if compressed_path.exists():
                            logger.info(f"FITS file already exists at {compressed_path.exists()}. Skipping.")
                            continue
                        try:
                            if mode == "proj":
                                yt_fits = FITSOffAxisProjection(
                                    ds,
                                    normal=axis,
                                    fields=(
                                        "gas",
                                        f"xray_photon_intensity_" f"{energies.emin}_{energies.emax}_keV",
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
                                        f"xray_photon_intensity_" f"{energies.emin}_{energies.emax}_keV",
                                    ),
                                    center=(sub["pos_x"], sub["pos_y"], sub["pos_z"]),
                                    width=w,
                                    image_res=r,
                                )

                            yt_fits.update_header(field="all", key="AXIS", value=f"{axis}")
                            yt_fits.update_header(field="all", key="WIDTH", value=w)
                            yt_fits.update_header(field="all", key="REDSHIFT", value=redshift)
                            yt_fits.update_header(field="all", key="EMIN", value=energies.emin)
                            yt_fits.update_header(field="all", key="EMAX", value=energies.emax)
                            yt_fits.writeto(f"{fits_path}")
                            compress_gzip(in_file_path=fits_path, out_file_path=compressed_path, remove_file=True)
                        except Exception as e:
                            if environment.fail_on_error:
                                logger.exception(f"Failed to process {cutout.resolve()} with error:\n" f"{e}")
                                raise
                            else:
                                logger.warning(f"Failed to process {cutout.resolve()} with error:\n" f"{e}")
                        fits.append(compressed_path)
    if environment.consume_data:
        cutout.unlink()
        logger.success(f"Deleted {cutout}.")
    end = datetime.now()
    logger.info(f"Processing of {cutout} took {end - start}.")
    return fits
