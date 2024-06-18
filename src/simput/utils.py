import os
import shutil
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from loguru import logger
from xspec import Model, Xset
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table

from src.sixte import commands


def get_spectrumfile(run_dir: Path, norm=0.01) -> Path:
    spectrum_file = run_dir / "spectrum.xcm"
    if not spectrum_file.exists():
        logger.info(f"Spectrum file at {spectrum_file.resolve()} does not exist. Will create a new one.")

        with open(os.devnull, "w") as f, redirect_stdout(f), redirect_stderr(f):
            Model("phabs*power", setPars={1: 0.04, 2: 2.0, 3: norm})
            Xset.save(f"{spectrum_file.resolve()}")

    return spectrum_file


def merge_simputs(simput_files: list[Path], output_file: Path) -> Path:
    # Combine the simput point sources
    if len(simput_files) == 1:
        file = simput_files[0]
        shutil.copy2(file, output_file)
    else:
        commands.simputmerge(infiles=simput_files, outfile=output_file, fetch_extension=True)

    return output_file

    
def add_columns_to_simput(simput_file: Path, additional_columns: dict[str, np.ndarray]) -> None:
    
    with fits.open(simput_file, mode='update', memmap=False) as hdul:
        # Identify the extension with source data
        source_hdu = None
        for hdu in hdul:
            if isinstance(hdu, fits.BinTableHDU) and 'RA' in hdu.columns.names and 'DEC' in hdu.columns.names:
                source_hdu = hdu
                break

        if source_hdu is None:
            raise ValueError("No suitable extension with RA and DEC columns found in the SIMPUT file.")

        # Convert additional columns to fits.Columns
        new_columns = []
        for col_name, col_data in additional_columns.items():
            col_format = 'D' if np.issubdtype(col_data.dtype, np.floating) else 'J'  # Use 'D' for floats, 'J' for integers
            new_col = fits.Column(name=col_name, format=col_format, array=np.array(col_data))  # Force load data
            new_columns.append(new_col)

        # Convert existing columns to fits.Columns
        existing_columns = []
        for col in source_hdu.columns:
            col_data = np.array(source_hdu.data[col.name])  # Force load data
            existing_col = fits.Column(name=col.name, format=col.format, array=col_data)
            existing_columns.append(existing_col)

        # Combine existing columns with new columns
        all_columns = fits.ColDefs(existing_columns + new_columns)

        # Create a new BinTableHDU with all columns
        new_hdu = fits.BinTableHDU.from_columns(all_columns)

        # Replace the old HDU with the new one
        hdul[hdul.index_of(source_hdu)] = new_hdu
        hdul.flush()
        
        
def make_deblending_img_settings(i, n_blended, img_settings):
    """Change the 'deblending' parameter in 'img_settings' to True if the index of the current image is smaller than the absolute number of images that shall contain blended sources.  

    Args:
        i (int): index of the current image
        n_blended (_type_): absolute number of images containing blended sources 

    Returns:
        _type_: 
    """                

    if i<n_blended:
        img_settings.update({"deblending":True})
        
    return img_settings
