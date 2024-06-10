from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal
from uuid import uuid4

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


def create_agn_sources(
    emin: float,
    emax: float,
    run_dir: Path,
    img_settings: dict,
    xspec_file: Path,
    offset: tuple[float, float] | str = (0.0, 0.0),
    center_point: tuple[float, float] = (0.0, 0.0),
):
    output_files = []
    
    
    #TODO: implement a vesion that works outside of the debugging mode 
    if img_settings["deblending_n_gen"] > 0:
        # Compute absolute number of images that should contain blended sources 
        abs_deblending_n_gen = int(img_settings["deblending_n_gen"]*img_settings["n_gen"])
           
    else:
        # contains_blended_sources = False 
        frac_blended_sources = 0 
        blended_offset_vals = None
        blended_idx = None
    
    #TODO: define proper path 
    # with open('deblending_stats.csv', mode='w', newline='') as file:
        
        # Create a CSV writer object
        # writer = csv.writer(file)

        # Write header row
        # writer.writerow(['idx','image_ID', 'contains_blended_sources', 'frac_blended_sources', 'fluxes', 'locations (ra, dec) [arcsec]', 'blended_idx', 'blended_fluxes1', 'blended_fluxes2', 'blended locations 1 (ra, dec) [arcsec]', 'blended_locations2 (ra, dec) [arcsec]', 'blended_offset_vals (ra, dec) [arcsec]'])
        # writer.writerow(['idx','image_ID', 'contains_blended_sources', 'frac_blended_sources', 'fluxes', 'ra [arcsec]', 'dec [arcsec]' ,'blended_idx0', 'blended_idx1', 'blended_fluxes0', 'blended_fluxes1', 'blended0 ra [arcsec]', 'blended1 ra [arcsec]', 'blended0 dec [arcsec]', 'blended1 dec [arcsec]', 'offset ra [arcsec]', 'offset dec [arcsec]'])
        
    for i_gen in range(img_settings["n_gen"]):
        # Use the current time as id, such that clashes don't happen
        unique_id = uuid4().int
        output_file_path = run_dir / f"agn_{unique_id}_p0_{emin}ev_p1_{emax}ev.simput"
        simput_files: list[Path] = []
        
        # Get the fluxes from the agn distribution
        fluxes = get_fluxes(img_settings["agn_counts_file"])
        num_fluxes = len(fluxes)
    
        # Compute the offsets 
        # The FOV is the same for EPN, EMOS1, and EMOS2
        fov = get_fov("epn")

        # Randomly position the point source within the fov
        rng = np.random.default_rng()
        if offset == "random":
            offset_vals = rng.uniform(low=-1.0 * fov / 2, high=fov / 2, size=(num_fluxes,2))
            
            
            if img_settings["put_source_in_center"]:
                # Force one agn to be at the center
                offset_vals[0] = np.array([0,0])
                fluxes [0] = 0.8*np.min(fluxes)  # Make it bright, so it is obvious which one is the one at the center
                num_fluxes-= 1                   # Do not count the center source in the number of fluxes

        
        
        # Making blended sources starting at 1 to not duplicate AGN in center of image
        if i_gen < abs_deblending_n_gen:
            
            # contains_blended_sources = True
            
            # Determine fraction of blended fluxes 
            frac_blended_sources = rng.normal(loc = img_settings["deblending_n_flux"], scale = 0.1)
            # Make sure that the value is in range [0,1]
            frac_blended_sources = np.clip(frac_blended_sources, 0, 0.5)
            
            # Compute absolute number of AGNs that should be blended 
            abs_deblending_n_flux = int(frac_blended_sources*(num_fluxes))
            
            # Determine indices of blended sources, starting from index 1 because index 0 is supposed to be at the center and should not be blended (for now)
            blended_idx = rng.choice(np.arange(int(img_settings["put_source_in_center"]), len(fluxes)), size = (abs_deblending_n_flux,2), replace = False)
            
            # Compute the offsets of the blended sources, divding by 3600 to go from arcseconds to degrees
            blended_dist = rng.uniform(low=img_settings["deblending_min_sep"]/3600, high=img_settings["deblending_max_sep"]/3600, size= (abs_deblending_n_flux, 2))
            # blended_dist = np.array([8/3600, 8/3600])
            blended_offset_vals = offset_vals[blended_idx[:,0]] + blended_dist
            
            # Change offset values of the blended sources that are moved 
            offset_vals[blended_idx[:,1]] = blended_offset_vals
            
        else:
            contains_blended_sources = False
            frac_blended_sources = 0
        
        location = (center_point[0] + offset_vals[:,0], center_point[1] + offset_vals[:,1])
        # loc_blend = np.array(location)[:,blended_idx[:,0]]
        # loc_blend2 = np.array(location)[:,blended_idx[:,1]]
        # flux_blend = fluxes[blended_idx[:,0]]
        # flux_blend2 = fluxes[blended_idx[:,1]]
        
        #TODO: implement more efficient way of logging 
        # Write data to CSV file
        
        # writer.writerow([i_gen, f"agn_{unique_id}_p0_{emin}ev_p1_{emax}ev", contains_blended_sources, frac_blended_sources, fluxes, np.transpose(np.array(location))*3600, blended_idx, flux_blend, flux_blend2, loc_blend*3600, loc_blend2*3600, blended_dist*3600])
        # row = [i_gen, f"agn_{unique_id}_p0_{emin}ev_p1_{emax}ev", contains_blended_sources, frac_blended_sources, fluxes, location[0]*3600, location[1]*3600, blended_idx[:,0], blended_idx[:,1], flux_blend, flux_blend2, loc_blend[0]*3600, loc_blend2[0]*3600, loc_blend[1]*3600, loc_blend2[1]*3600, blended_dist[:,0]*3600, blended_dist[:,1]*3600]
        
        # Create array that indicates the index of the partner blended sources for each source. Sources that are not blended are mapped onto
        deblending_indices = np.arange(len(fluxes))
        deblending_indices[blended_idx[:,0]] = blended_idx[:,1]
        deblending_indices[blended_idx[:,1]] = blended_idx[:,0]
        
        deblend_columns = {
            'DEBLENDING INDICES': deblending_indices,
            }
        
        # serialized_row = [json.dumps(cell.tolist()) if isinstance(cell, np.ndarray) else cell for cell in row]
        # writer.writerow(serialized_row)
        
        for i, flux in enumerate(fluxes):
            
            logger.info(f"Creating AGN with flux={flux}")
            output_file = run_dir / f"ps_{unique_id}_{i}.simput"
            output_file = simput_ps(
                emin=emin,
                emax=emax,
                output_file=output_file,
                location = (location[0][i],location[1][i]),
                src_flux=flux,
                offset = offset_vals[i],
                xspec_file=xspec_file,
            )
            simput_files.append(output_file)
        output_file = merge_simputs(simput_files=simput_files, output_file=output_file_path)
        
        # Add additional information for deblending
        # add_columns_to_simput(output_file, deblend_columns)
        
        output_files.append(output_file)
    

        for file in simput_files:
            file.unlink(missing_ok=True)

    return output_files


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

        if mode == "agn":
            file_names = create_agn_sources(
                emin=emin,
                emax=emax,
                run_dir=run_dir,
                img_settings=img_settings,
                xspec_file=spectrum_file,
                offset = "random"
            )

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
