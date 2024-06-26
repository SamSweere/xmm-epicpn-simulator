
import tarfile
import io
from astropy.io import fits
import os 
import numpy as np 
import shutil


# Function to find all .fits files in a specific subfolder
def find_fits_files(directory):
    fits_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.fits.gz'):
                fits_files.append(os.path.join(root, file))
    return fits_files


def make_ds9_region_file(hdul):
    
    img_hdu = hdul[0]
    agn_hdu = hdul[1]
    
    # idx = agn_hdu["SRC_ID"]
    # image_ID = agn_hdu["image_ID"]
    
    ra_locs = agn_hdu.data["RA"]
    dec_locs = agn_hdu.data["DEC"]
    
    deblending_indices = agn_hdu.data["DEBLENDING INDICES"]
    
    center_ra = ra_locs[0:1]
    center_dec = dec_locs[0:1]
    
    # Find the single sources 
    single_idx = np.where(np.arange(len(ra_locs))== np.array(deblending_indices))[0][1:]
    blended_idx = np.where(np.arange(len(ra_locs))!= np.array(deblending_indices))[0]
    
    types = [
    {
        "x_positions": ra_locs[single_idx],  # X coordinates for the first type
        "y_positions": dec_locs[single_idx],  # Y coordinates for the first type
        "color": "black",  # Color for the first type
        "text_template": "non-blended sources"  # Text template for the first type
    },
    {
        "x_positions": ra_locs[blended_idx],  # X coordinates for the second type
        "y_positions": dec_locs[blended_idx],  # Y coordinates for the second type
        "color": "blue",  # Color for the second type
        "text_template": "blended sources"  # Text template for the second type
    },
    {
        "x_positions": center_ra,  # X coordinates for the second type
        "y_positions": center_dec,  # Y coordinates for the second type
        "color": "red",  # Color for the second type
        "text_template": "center_source"  # Text template for the second type
    },
    
    ]
        
    # Common parameter for the circles
    radius = '20"'  # Radius of the circles

    # Specify the file name
    simput_name = img_hdu.header["SIMPUT"]
    #Split the filename twice to remove both extensions
    root, ext1 = os.path.splitext(simput_name)
    root, ext2 = os.path.splitext(root)
    file_name = f'data/xmm_sim_dataset/epn/thin/agn/{root}.reg'

    # Define the content of the DS9 regions file
    header = "fk5\n"

    regions = []
    for type_index, type_info in enumerate(types):
        x_positions = type_info["x_positions"]
        y_positions = type_info["y_positions"]
        color = type_info["color"]
        text = type_info["text_template"]
        
        for i, (x, y) in enumerate(zip(x_positions, y_positions)):
           
            region = f"circle({x}, {y}, {radius}) #color={color} text={{{text}}}"
            region = f"circle({x}, {y}, {radius}) #color={color}"
            regions.append(region)

    # Write the content to the .reg file
    with open(file_name, 'w') as reg_file:
        # Write the header
        reg_file.write(header)
        # Write each region
        for region in regions:
            reg_file.write(region + "\n")

    print(f"{file_name} has been created successfully.")
    
    
    with tarfile.open(tar_gz_path, 'w:gz') as tar:
        tar.add(extracted_folder+'/agn', arcname = '')
        
    # Delete the uncompressed folder
    shutil.rmtree(extracted_folder)
    



# Path to the tar.gz file
tar_gz_path = 'data/xmm_sim_dataset/epn/thin/agn.tar.gz'
subfolder = 'agn/10ks/1x/'

# Extract the tar.gz file
extracted_folder = 'data/xmm_sim_dataset/epn/thin'
with tarfile.open(tar_gz_path, 'r:gz') as tar:
    tar.extractall(path=extracted_folder)

# Define the specific subfolder within the extracted files
specific_subfolder = os.path.join(extracted_folder, subfolder)

# Get all .fits files in the specific subfolder
fits_files = find_fits_files(specific_subfolder)

# Perform operations on each .fits file
for fits_file in fits_files:
    with fits.open(fits_file) as hdul:
        # Perform operations on the FITS file
        print(f"Processing file: {fits_file}")
        
        test = hdul[1].columns.names
        test0 = hdul[0].header
        
        # Make the region file
        make_ds9_region_file(hdul)
        
