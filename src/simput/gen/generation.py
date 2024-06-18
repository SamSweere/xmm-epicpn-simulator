from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal
from uuid import uuid4

import numpy as np
from loguru import logger

from src.simput.agn import get_fluxes
from src.simput.gen.background import background
from src.simput.gen.image import simput_image
from src.simput.gen.pointsource import simput_ps
from src.simput.utils import merge_simputs
from src.xmm_utils.file_utils import compress_gzip
from src.xmm.utils import get_fov
import numpy as np
from src.simput.utils import add_columns_to_simput

import csv
import json




def create_background(
    instrument_name: Literal["epn", "emos1", "emos2"],
    emin: float,
    emax: float,
    run_dir: Path,
    spectrum_file: Path,
) -> list[Path]:
    output_files = [
        background(
            run_dir=run_dir,
            spectrum_file=spectrum_file,
            instrument_name=instrument_name,
            emin=emin,
            emax=emax,
        )
    ]

    return output_files


def create_agn_simput(
    agn_counts_file: Path,
    emin: float,
    emax: float,
    run_dir: Path,
    img_settings: dict,
    xspec_file: Path,
    offset: tuple[float, float] | str = (0.0, 0.0),
)-> None:
    fov = img_settings["fov"]
    instruments = img_settings["instruments"]
    center_points = img_settings["center_points"]
    output_dirs = img_settings["output_dirs"]
    simput_agn_settings = img_settings["simput_agn_settings"]
    deblending = img_settings["deblending"]
    
    if isinstance(offset, str) and offset != "random":
        raise ValueError(f'Value of offset is unknown string "{offset}"!')

    rng = np.random.default_rng()

    # Use the current time as id, such that clashes don't happen
    unique_id = uuid4().int
    final_name = f"agn_{unique_id}_p0_{emin}ev_p1_{emax}ev.simput.gz"
    simput_files: list[Path] = []
    
    # Get the fluxes from the agn distribution
    fluxes = get_fluxes(agn_counts_file)
    num_fluxes = fluxes.shape[0]

    if isinstance(offset, str) and offset == "random":
        offset_vals = rng.uniform(low=-fov / 2.0, high=fov / 2.0, size=(num_fluxes, 2))
        
    if simput_agn_settings.put_source_in_center:
        # Force one agn to be at the center
        offset_vals[0] = np.array([0,0])
        fluxes [0] = 0.8*np.min(fluxes)  # Make it less bright than other agns 
        num_fluxes-= 1                   # Do not count the center source in the number of fluxes

    
    # Making blended sources starting at 1 to not duplicate AGN in center of image
    if deblending:
        
        # Determine fraction of blended fluxes 
        frac_blended_sources = rng.normal(loc = simput_agn_settings.deblending_n_flux, scale = 0.1)
        
        # Make sure that the value is in range [0,1]
        frac_blended_sources = np.clip(frac_blended_sources, 0, 0.5)
        
        # Compute absolute number of AGNs that should be blended 
        abs_deblending_n_flux = int(frac_blended_sources*(num_fluxes))
        
        # Determine indices of blended sources, starting from index 1 because index 0 is supposed to be at the center and should not be blended (for now)
        blended_idx = rng.choice(np.arange(int(simput_agn_settings.put_source_in_center), len(fluxes)), size = (abs_deblending_n_flux,2), replace = False)
        
        # Compute the offsets of the blended sources, divding by 3600 to go from arcseconds to degrees
        blended_dist = rng.uniform(low=simput_agn_settings.deblending_min_sep/3600, high=simput_agn_settings.deblending_max_sep/3600, size= (abs_deblending_n_flux, 2))
        # blended_dist = np.array([8/3600, 8/3600])
        blended_offset_vals = offset_vals[blended_idx[:,0]] + blended_dist
        
        # Change offset values of the blended sources that are moved 
        offset_vals[blended_idx[:,1]] = blended_offset_vals
        
    else:
        contains_blended_sources = False
        frac_blended_sources = 0
    
    location = (center_points[0][0] + offset_vals[:,0], center_points[0][1] + offset_vals[:,1])
   
    # Create array that indicates the index of the partner blended sources for each source. Sources that are not blended are mapped onto
    deblending_indices = np.arange(len(fluxes))
    deblending_indices[blended_idx[:,0]] = blended_idx[:,1]
    deblending_indices[blended_idx[:,1]] = blended_idx[:,0]
    
    deblend_columns = {
        'DEBLENDING INDICES': deblending_indices,
        }
    
    for instrument, output_dir in zip(instruments, output_dirs, strict=False):
        logger.info(f"Creating AGN SIMPUTs for {instrument} with id {unique_id}.")
        
        for i, flux in enumerate(fluxes):
            logger.debug(f"Creating AGN with flux={flux}")
            output_file = run_dir / f"ps_{unique_id}_{i}.simput"
            output_file = simput_ps(
                emin=emin,
                emax=emax,
                output_file=output_file,
                location = (location[0][i],location[1][i]),
                src_flux=flux,
                xspec_file=xspec_file,
            )
            simput_files.append(output_file)
        merged = merge_simputs(simput_files=simput_files, output_file=run_dir / f"merged_{unique_id}.simput")
        
        # Add additional information for deblending
        add_columns_to_simput(output_file, deblend_columns)
        
        compress_gzip(in_file_path=merged, out_file_path=output_dir / final_name, remove_file=True)
    

    for file in simput_files:
        file.unlink(missing_ok=True)


def simput_generate(
    emin: float,
    emax: float,
    mode: str,
    img_settings: dict,
    tmp_dir: Path,
    output_dir: Path,
    spectrum_file: Path,
) -> None:
    with TemporaryDirectory(dir=tmp_dir) as temp:
        run_dir = Path(temp)

        file_names = []

        if mode == "img":
            file_names = simput_image(
                emin=emin,
                emax=emax,
                run_dir=run_dir,
                img_settings=img_settings,
                xspec_file=spectrum_file,
            )
            img_path_in = img_settings["img_path"]
            tng_name = img_path_in.parts[-3]
            snapshot_num = img_path_in.parts[-2]
            output_dir = output_dir / tng_name / snapshot_num
            output_dir.mkdir(parents=True, exist_ok=True)

        if mode == "bkg":
            file_names = create_background(
                instrument_name=img_settings["instrument_name"],
                emin=emin,
                emax=emax,
                run_dir=run_dir,
                spectrum_file=spectrum_file,
            )

        for file_name in file_names:
            # Compress the simput file and move it to the correct output dir
            compressed_file = output_dir / f"{file_name.name}.gz"
            if compressed_file.exists():
                logger.warning(f"SIMPUT file {compressed_file.resolve()} already exists, skipping.")
            else:
                compress_gzip(in_file_path=file_name, out_file_path=compressed_file)
            file_name.unlink(missing_ok=True)
