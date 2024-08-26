
import astropy 
import numpy as np
import glob
from astropy.io import fits
import os 

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.optimize import curve_fit
from matplotlib.colors import LogNorm


def remove_duplicates(pairs):
    """
    Remove duplicate pairs from an array where pairs like (x, y) and (y, x) are considered the same.
    
    Parameters:
    pairs (np.ndarray): An array of shape (n, 2) containing the pairs.
    
    Returns:
    np.ndarray: An array containing unique pairs.
    """
    # Step 1: Convert each pair to a sorted tuple
    sorted_pairs = [tuple(sorted(pair)) for pair in pairs]
    
    # Step 2: Use a set to store unique pairs
    unique_pairs_set = set(sorted_pairs)
    
    # Step 3: Convert the set back to a NumPy array
    unique_pairs_array = np.array(list(unique_pairs_set))
    
    return unique_pairs_array


def reverse_access(array):
    """
    Create an array where the values are indices and the indices are values of the original array.
    
    Parameters:
    array (np.ndarray): The original array of values.
    
    Returns:
    np.ndarray: The array with reversed access.
    """
    # Initialize the new array with -1 or a suitable default value
    # The size of the new array should be at least max(array) + 1
    reversed_array = np.full(np.max(array) + 1, -1)
    
    # Iterate through the original array and set the values in the reversed array
    for index, value in enumerate(array):
        reversed_array[value] = index
    
    return reversed_array


def find_bounding_boxes(array, buffer = 5):
    """
    Find bounding boxes for each pair in the given array.
    
    Parameters:
    array (np.ndarray): An array of shape (2, 11, 2) where the dimensions refer to:
                        2 (coordinates: x and y), 11 (number of point pairs), 2 (point pairs)
    
    Returns:
    np.ndarray: An array of shape (11, 4) where each row contains the bounding box [xmin, ymin, xmax, ymax]
    """
   

    # Extract the x and y coordinates for the pair
    x_coords = array[0,:]
    y_coords = array[1,:]
    
    # Compute the bounding box
    xmin = np.min(x_coords) - buffer
    xmax = np.max(x_coords) + buffer
    ymin = np.min(y_coords) - buffer
    ymax = np.max(y_coords) + buffer
    
    # Store the bounding box
    bounding_box = [xmin, ymin, xmax, ymax]
    
    return bounding_box

# TODO: remove one of the sigmas 
def gaussian_2d(xy, amp, x0, y0, sigma_x, sigma_y, theta, offset):
    """
    2D Gaussian function.
    
    Parameters:
    xy (tuple): x and y coordinates.
    amp (float): Amplitude of the Gaussian.
    x0 (float): x-coordinate of the center.
    y0 (float): y-coordinate of the center.
    sigma_x (float): Standard deviation in the x direction.
    sigma_y (float): Standard deviation in the y direction.
    theta (float): Rotation angle in radians.
    offset (float): Offset of the Gaussian.
    
    Returns:
    np.ndarray: Gaussian values at the given coordinates.
    """
    x, y = xy
    a = np.cos(theta)**2 / (2 * sigma_x**2) + np.sin(theta)**2 / (2 * sigma_y**2)
    b = -np.sin(2 * theta) / (4 * sigma_x**2) + np.sin(2 * theta) / (4 * sigma_y**2)
    c = np.sin(theta)**2 / (2 * sigma_x**2) + np.cos(theta)**2 / (2 * sigma_y**2)
    g = offset + amp * np.exp(- (a * (x - x0)**2 + 2 * b * (x - x0) * (y - y0) + c * (y - y0)**2))
    return g

def fit_gaussian(image, point_sources):
    """
    Fit 2D Gaussians to the given point sources in the image.
    
    Parameters:
    image (np.ndarray): The input image array.
    point_sources (list of tuples): Known positions of the point sources.
    
    Returns:
    list of tuples: Fitted parameters for each Gaussian.
    """
    def fit_single_gaussian(image, point_source):
        """
        Fit a single 2D Gaussian to the point source.
        
        Parameters:
        image (np.ndarray): The cropped image.
        point_source (tuple): Known position of the point source.
        
        Returns:
        tuple: Fitted parameters of the Gaussian.
        """
        x0, y0 = point_source
        x = np.arange(image.shape[1])
        y = np.arange(image.shape[0])
        x, y = np.meshgrid(x, y)
        initial_guess = (image[y0, x0], x0, y0, 1, 1, 0, 0)
        
        def gaussian_2d_flat(xy, amp, x0, y0, sigma_x, sigma_y, theta, offset):
            x, y = xy
            return gaussian_2d((x, y), amp, x0, y0, sigma_x, sigma_y, theta, offset)
        
        #TODO: check if you can add constraints to fitting parameters 
        popt, _ = curve_fit(gaussian_2d_flat, (x.ravel(), y.ravel()), image.ravel(), p0=initial_guess)
       
        return popt
    
    params = []
    for point_source in point_sources:
        params.append(fit_single_gaussian(image, point_source))
    
    return params

   
def plot_results(image, point_sources, params, SRC_IDs, pixel_coords):
    """
    Plot the cropped image and fitted Gaussians using logarithmic scale.
    
    Parameters:
    image (np.ndarray): The cropped image.
    point_sources (list of tuples): Known positions of the point sources.
    params (list of tuples): Fitted parameters of the Gaussians.
    """
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))

    # Normalize the original image
    image_min = np.min(image[image > 0])  # Avoid log of zero
    image_max = np.max(image)
    
    # Plot original image with logarithmic scale
    im1 = axs[0].imshow(image, cmap='viridis', norm=LogNorm(vmin=image_min, vmax=image_max))
    axs[0].set_title('Original Image')
    axs[0].set_xlabel('Pixel X')
    axs[0].set_ylabel('Pixel Y')
    axs[0].axis('on')
    
    # Generate a higher resolution grid for plotting the Gaussian fits
    x = np.linspace(0, image.shape[1], 5*image.shape[1])
    y = np.linspace(0, image.shape[0], 5*image.shape[0])
    x, y = np.meshgrid(x, y)
    combined_gaussians = np.zeros_like(x, dtype=float)
    
    # Plot each Gaussian fit
    for i, (param, point_source) in enumerate(zip(params, point_sources)):
        amp, x0, y0, sigma_x, sigma_y, theta, offset = param
        fitted_gaussian = gaussian_2d((x, y), amp, x0, y0, sigma_x, sigma_y, theta, offset)
        combined_gaussians += fitted_gaussian
    
    # Normalize the combined Gaussians
    gaussian_min = np.min(combined_gaussians[combined_gaussians > 0])  # Avoid log of zero
    gaussian_max = np.max(combined_gaussians)
    
    # Display the combined Gaussians with logarithmic scale
    im2 = axs[1].imshow(combined_gaussians, extent=[0, image.shape[1], 0, image.shape[0]], origin='lower', cmap='viridis', norm=LogNorm(vmin=gaussian_min, vmax=gaussian_max), interpolation='bilinear')
    axs[1].set_title('Fitted 2D Gaussians')
    axs[1].set_xlabel('Pixel X')
    axs[1].set_ylabel('Pixel Y')
    axs[1].axis('on')
    
    fig.savefig(f'data/gaussian_fit_images/gaussian_fit_ID_{SRC_IDs}_coords_{pixel_coords}.pdf')

#Define the basepath
basepath = 'data/xmm_sim_dataset/epn/thin/agn'
# Get the names of all the fits files
lr_files  = glob.glob(basepath +'/20ks/1x/agn_*_p_0-4*') # --> Try with join
hr_files = glob.glob(basepath + '/100ks/2x/agn_*')


for hr_file in hr_files:
    
    # Load the fits file: 
    with fits.open(hr_file) as hdul:
    
        primary_hdu = hdul[0]
        source_catalogue = hdul[1]

        # Access the data
        img = primary_hdu.data
        img_header = primary_hdu.header
        source_catalogue_data = source_catalogue.data
        
        # Plot the whole image 
        fig, ax = plt.subplots()
        img[img <= 0] = np.min(img[img>0])
        ax.imshow(img, norm = LogNorm(vmin = np.min(img[img>0]), vmax = np.max(img)))
        fig.savefig('data/gaussian_fit_images/whole_image.pdf')
        
        SRC_IDs = source_catalogue_data['SRC_ID'] -1
        # Reverse the access to the array such that we can find out the index of a specific source ID
        idx_SRC_ID = reverse_access(SRC_IDs)
        
        deblending_indices = source_catalogue_data['DEBLENDING INDICES']
        
        RA = source_catalogue_data['RA']
        DEC = source_catalogue_data['DEC']
        
        # Find the blended sources
        blended_idx = np.where(SRC_IDs!= deblending_indices)
        
        ID_pairs = np.stack((SRC_IDs[blended_idx], deblending_indices[blended_idx]), axis = 1)
        ID_unique_pairs = remove_duplicates(ID_pairs)
        idx_unique_pairs = idx_SRC_ID[ID_unique_pairs]
        
        # Find the positions of the pairs 
        blended_RA = RA[idx_unique_pairs]
        blended_DEC = DEC[idx_unique_pairs]
        
        # Convert the positions to pixel coordinates 
        
        wcs = WCS(img_header)
    
        # Convert AGN fk5 coordinates to pixel coordinates
        agn_coords = SkyCoord(ra=blended_RA*u.degree, dec=blended_DEC*u.degree, frame='fk5')
        pixel_coords = np.array(agn_coords.to_pixel(wcs)).astype(int)
        
        # Loop through all the blended pairs
        for i in range(pixel_coords.shape[1]):
            
            # Find the bounding box for the current pair
            coords = pixel_coords[:, i, :]
            bounding_box = find_bounding_boxes(coords, buffer = 8)
            
            # Compute pixel coorindates within new, cropped frame
            x1 = coords[0,0] - bounding_box[0]
            y1 = coords[1,0] - bounding_box[1]
            x2 = coords[0,1] - bounding_box[0]
            y2 = coords[1,1] - bounding_box[1]
            point_sources = [(x1, y1), (x2, y2)]
       
            # Crop the current pair from the image 
            cropped_img = img[bounding_box[1]:bounding_box[3], bounding_box[0]:bounding_box[2]]
            
            if not np.all(img == 0):
                # Fit the Gaussians
                params = fit_gaussian(cropped_img, point_sources)

                # Plot the results
                plot_results(cropped_img, point_sources, params, SRC_IDs=ID_unique_pairs[i], pixel_coords = coords[:,0])
           
    
        test = 5
        
        # Determine a window to cut 
        
    

# Make the environment 
# Install astropy 

# Load a fits file 
# function that gets crop from position of AGNs --> which AGN? and do I get all the crops?
# So I am looping through all AGNs and all images?


# some function 
# input: crop of region 
# fits a 2D-gaussian to the two sources --> Do we really need the Gaussian? Can't we just apply it directly to the data?
# Find some sort of local minima between the two 
# Compute the measure 

# Save the results of the measure somehow 


# Do this for all low resolution images 
# Do this for all high resolution images 

# Make a plot showing the measure as a function of the distance between the two sources also account for position (but what do I use for the position? the minimum inbetween them?) of the AGNs?
# Could also make a plot showing the distribution of the distances 


