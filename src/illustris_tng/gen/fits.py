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
        width: List[int] = None,
        resolution: List[int] = None,
        redshift=0.05,
        overwrite=False
):
    if width is None:
        width = [1000]
    if resolution is None:
        resolution = [2048]
    try:
        logger.info(f"Processing cutout: {cutout.resolve()}. Arguments: mode={mode_dict}; emin={emin}; emax={emax}; "
                    f"width={width}; resolution={resolution}; redshift={redshift}; overwrite={overwrite}; "
                    f"output_dir={output_dir}; cutout_datafolder={cutout.parent}")
        start = datetime.now()
        ds = yt.load(f"{cutout.resolve()}", default_species_fields="ionized")

        yt.add_xray_emissivity_field(ds, emin, emax, redshift, data_dir=cloudy_emissivity_root)

        for mode, normals in mode_dict.items():
            for normal in normals:
                for w in width:
                    for r in resolution:
                        normal_str = f"{normal[0]}_{normal[1]}_{normal[2]}" if mode == "proj" else normal
                        filename = f"{cutout.name.replace('.hdf5', '')}_m_{mode}_r_{r}_w_{w}_n_{normal_str}"

                        fits_filename = filename + '.fits'
                        fits_path = output_dir / fits_filename

                        if fits_path.exists() and not overwrite:
                            raise FileExistsError(f"{fits_path} already exists and `overwrite` is False")
                        if mode == "proj":
                            yt_fits = yt.FITSOffAxisProjection(ds, normal=normal,
                                                               fields=('gas', 'xray_photon_intensity_0.5_2.0_keV'),
                                                               center=(sub['pos_x'], sub['pos_y'], sub['pos_z']),
                                                               width=(w, "code_length"), image_res=r)
                        else:
                            yt_fits = yt.FITSSlice(ds, axis=normal, fields=('gas', 'xray_photon_intensity_0.5_2.0_keV'),
                                                   center=(sub['pos_x'], sub['pos_y'], sub['pos_z']),
                                                   width=(w, "code_length"), image_res=r)

                        yt_fits.update_header(field="all", key="NORMAL", value=str(normal))
                        yt_fits.update_header(field="all", key="WIDTH", value=w)
                        yt_fits.update_header(field="all", key="REDSHIFT", value=redshift)
                        yt_fits.update_header(field="all", key="EMIN", value=emin)
                        yt_fits.update_header(field="all", key="EMAX", value=emax)

                        yt_fits.writeto(str(fits_path.resolve()), overwrite=overwrite)
        end = datetime.now()
        logger.info(f"Processing of {cutout} done!")
        logger.info(f"Duration: {end - start}")

    except Exception as e:
        logger.exception(f" ERROR, failed to process {cutout.name} with errro: {e}")
