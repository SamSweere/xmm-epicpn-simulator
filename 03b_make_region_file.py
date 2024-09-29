import tarfile
import os
import numpy as np
from astropy.io import fits
from argparse import ArgumentParser


def find_fits_files(directory):
    """
    Find all .fits.gz files in the specified directory and its subdirectories.

    Parameters:
    - directory: str, the path to the directory to search for .fits.gz files.

    Returns:
    - List[str]: A list of paths to the .fits.gz files found in the directory.
    """
    fits_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.fits.gz'):
                fits_files.append(os.path.join(root, file))
    return fits_files


def make_ds9_region_file(hdul, tar_gz_path):
    """
    Create a DS9 region file from the provided HDU list and package it into a tar.gz file.

    Parameters:
    - hdul: HDU list containing FITS data.
    - tar_gz_path: Path to the original tar.gz file for naming the output.
    """
    img_hdu = hdul[0]
    agn_hdu = hdul[1]

    # Extract relevant data
    ra_locs = agn_hdu.data["RA"]
    dec_locs = agn_hdu.data["DEC"]
    deblending_indices = agn_hdu.data["DEBLENDING INDICES"]
    SRC_IDs = agn_hdu.data["SRC_ID"] - 1  # Adjust for zero-based indexing

    center_ra = ra_locs[0:1]
    center_dec = dec_locs[0:1]

    # Identify single and blended sources
    single_idx = np.where(SRC_IDs == np.array(deblending_indices))[0][1:]
    blended_idx = np.where(SRC_IDs != np.array(deblending_indices))[0]

    # Prepare region types
    types = [
        {
            "x_positions": ra_locs[blended_idx],
            "y_positions": dec_locs[blended_idx],
            "color": "blue",
            "text_template": "blended sources"
        },
        {
            "x_positions": center_ra,
            "y_positions": center_dec,
            "color": "red",
            "text_template": "center_source"
        },
        {
            "x_positions": ra_locs[single_idx],
            "y_positions": dec_locs[single_idx],
            "color": "green",
            "text_template": "non-blended sources"
        },
    ]

    radius = '20"'  # Radius of the circles
    simput_name = img_hdu.header["SIMPUT"]
    root, ext1 = os.path.splitext(simput_name)
    root, ext2 = os.path.splitext(root)
    file_name = f'data/xmm_sim_dataset/epn/thin/agn/{root}.reg'

    # Define the content of the DS9 regions file
    header = "fk5\n"
    regions = []
    
    for type_info in types:
        x_positions = type_info["x_positions"]
        y_positions = type_info["y_positions"]
        color = type_info["color"]
        for x, y in zip(x_positions, y_positions):
            region = f"circle({x}, {y}, {radius}) #color={color}"
            regions.append(region)

    # Write the content to the .reg file
    with open(file_name, 'w') as reg_file:
        reg_file.write(header)
        for region in regions:
            reg_file.write(region + "\n")

    print(f"{file_name} has been created successfully.")

    # Create a new tar.gz file for the regions
    new_tar_gz_path = tar_gz_path.replace('agn', 'agn_region')
    with tarfile.open(new_tar_gz_path, 'w:gz') as tar:
        tar.add('data/xmm_sim_dataset/epn/thin/agn', arcname='')

    print(f"The zipped agn_region file can be found at: {new_tar_gz_path}")

def main(tar_gz_path, subfolder, extracted_folder):
    # Extract the tar.gz file
    with tarfile.open(tar_gz_path, 'r:gz') as tar:
        tar.extractall(path=extracted_folder)

    # Define the specific subfolder within the extracted files
    specific_subfolder = os.path.join(extracted_folder, subfolder)

    # Get all .fits files in the specific subfolder
    fits_files = find_fits_files(specific_subfolder)

    # Process each .fits file
    for fits_file in fits_files:
        with fits.open(fits_file) as hdul:
            print(f"Processing file: {fits_file}")
            make_ds9_region_file(hdul, tar_gz_path)

if __name__ == "__main__":
    parser = ArgumentParser(description="Create region files to mark blended AGN locations in DS9.")
    parser.add_argument(
        "-t", "--tar_gz_path",
        type=str,
        default='data/xmm_sim_dataset/epn/thin/agn.tar.gz',
        help="Path to the zipped tar.gz folder containing the fits files."
    )
    parser.add_argument(
        "-s", "--subfolder",
        type=str,
        default='agn/100ks/2x/',
        help="Subfolder for which region files shall be created. (The region files are the same for all exposures and resolutions, so it shouldnt really matter.)"
    )
    parser.add_argument(
        "-e", "--extracted_folder",
        type=str,
        default='data/xmm_sim_dataset/epn/thin',
        help="Folder where the tar.gz file will be extracted."
    )

    args = parser.parse_args()
    main(args.tar_gz_path, args.subfolder, args.extracted_folder)
