# This code generates the background noise (sky + instrument + particle) based on a background spectrum
import os

from astropy.io import fits
import numpy as np

from utils.external_run import run_headas_command


def create_exposure_fits(run_dir):
    # Generate the background image used for the simput
    n = np.ones((402, 402))
    hdu = fits.PrimaryHDU(n)

    hdul = fits.HDUList([hdu])

    header = hdul['PRIMARY'].header
    header['MTYPE1'] = "EQPOS"
    header['MFORM1'] = "RA,DEC"
    header['CTYPE1'] = "RA---TAN"
    header['CTYPE2'] = "DEC--TAN"

    header['CRPIX1'] = 186  # For some reason I need to offset it in the other direction#226 #n.shape[0] / 2  # 0.0
    header['CRPIX2'] = 219  # n.shape[1] / 2  # 0.0
    header['CRVAL1'] = 0.0
    header['CRVAL2'] = 0.0

    header['CUNIT1'] = "deg"
    header['CUNIT2'] = "deg"

    # cdelt give the pixel sizes in degrees
    # cdelt from XMM_MISCDATA_0022.CCF PLATE_SCALE_X, the unit is in arsec, arsec to degree by deciding it by 3600
    cdelt = (4.12838 / 3600)

    header['CDELT1'] = cdelt
    header['CDELT2'] = cdelt
    header['comment'] = "This fits image has all pixel value as 1 and has a similar resolution as xmm"

    out_file = os.path.join(run_dir, 'exposure_map.fits')
    hdul.writeto(out_file, overwrite=True)

    return out_file


def generate_ascii_spectrum(run_dir, spectrum_file, verbose):
    # Open the background spectrum file (sky + instrument + particle)
    hdu = fits.open(spectrum_file)

    bin_factor = hdu['SPECTRUM'].header['SPECDELT']
    channels = hdu['SPECTRUM'].data['CHANNEL']
    energies = channels * bin_factor / 1000

    # Calculate the rate based on the counts
    # counts = hdu['SPECTRUM'].data['COUNTS'].astype(np.float32)
    # rates = counts / float(hdu['SPECTRUM'].header['EXPOSURE'])

    pixel_size = 150e-4  # cm #based of: https://www.cosmos.esa.int/web/xmm-newton/boundaries-pn
    resolution = 402
    surface = (pixel_size * resolution) ** 2  # cm**2
    # cgi_rates = rates / surface  # photon/s/cm**2/keV

    out_file = os.path.join(run_dir, "exposure_spectrum.txt")
    with open(out_file, "w") as f:
        for i in range(len(energies)):
            rate = 1.0e-1  # *1.0e-4 #cgi_rates[i]
            energy_bin = energies[i]
            f.write(f"{energy_bin} {rate}")
            f.write("\n")

            if verbose:
                print(f"{energy_bin} keV: {rate} photon/s/cm**2/keV")

    if verbose:
        print(f"Ascii spectrum generated and saved to: {out_file}")

    return out_file


def generate_exposure_map_simput(run_dir, spectrum_file, emin, emax, verbose):
    # Create the background fits
    background_file = create_exposure_fits(run_dir)

    # Create the ascii spectrum file
    ascii_spectrum_file = generate_ascii_spectrum(run_dir, spectrum_file, verbose)

    simput_file_name = "exposure_map.simput"
    outfile_path = os.path.join(run_dir, simput_file_name)

    # Run the simputgen command
    simputgen_command = f"simputfile Simput={outfile_path} RA=0.0 DEC=0.0 Emin={emin} Emax={emax} history=True " \
                        f"clobber=True ImageFile={background_file} ASCIIFile={ascii_spectrum_file}"
    run_headas_command(simputgen_command, verbose=verbose)

    if verbose:
        print("Exposure generation complete")

    return simput_file_name
    # # Delete the background fits and ascii spectrum file
    # os.remove(background_file)
    # os.remove(ascii_spectrum_file)
