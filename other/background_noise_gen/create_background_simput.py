# This code generates the background noise (sky + instrument + particle) based on a background spectrum

from astropy.io import fits
import numpy as np

from utils import run_headas_command


# Generate the background image used for the simput
n = np.ones((402, 402))
hdu = fits.PrimaryHDU(n)

hdul = fits.HDUList([hdu])

header = hdul['PRIMARY'].header
header['MTYPE1'] = "EQPOS"
header['MFORM1'] = "RA,DEC"
header['CTYPE1'] = "RA---TAN"
header['CTYPE2'] = "DEC--TAN"

header['CRPIX1'] = n.shape[0]/2#0.0
header['CRPIX2'] = n.shape[1]/2#0.0
header['CRVAL1'] = 0.0
header['CRVAL2'] = 0.0

header['CUNIT1'] = "deg"
header['CUNIT2'] = "deg"

# cdelt give the pixel sizes in degrees
# cdelt from XMM_MISCDATA_0022.CCF PLATE_SCALE_X, the unit is in arsec, arsec to degree by deciding it by 3600
cdelt = (4.12838/3600)

header['CDELT1'] = cdelt
header['CDELT2'] = cdelt
header['comment'] = "This fits image has all pixel value as 1 and has a similar resolution as xmm"

hdul.writeto('const_background.fits', overwrite=True)




# Open the background spectrum file (sky + instrument + particle)
hdu = fits.open('pntffg_spectrum.ds')

bin_factor = hdu['SPECTRUM'].header['SPECDELT']
channels = hdu['SPECTRUM'].data['CHANNEL']
energies = channels*bin_factor/1000

# Calculate the rate based on the counts
counts = hdu['SPECTRUM'].data['COUNTS'].astype(np.float32)
rates = counts / float(hdu['SPECTRUM'].header['EXPOSURE'])

pixel_size = 150e-4 #cm #based of: https://www.cosmos.esa.int/web/xmm-newton/boundaries-pn
resolution = 402
surface = (pixel_size*resolution)**2 #cm**2
cgi_rates = rates/surface #photon/s/cm**2/keV

with open("background_spectrum.txt", "w") as f:
    for i in range(len(cgi_rates)):
        rate = cgi_rates[i]#*1.0e-4 #cgi_rates[i]
        energy_bin = energies[i]
        f.write(f"{energy_bin} {rate}")
        f.write("\n")

        print(f"{energy_bin} keV: {rate} photon/s/cm**2/keV")

simputgen_command = "simputfile Simput='background.simput' RA=0.0 DEC=0.0 Emin=0.15 Emax=15.0 history=True clobber=True ImageFile='const_background.fits' ASCIIFile='background_spectrum.txt'"
run_headas_command(simputgen_command)