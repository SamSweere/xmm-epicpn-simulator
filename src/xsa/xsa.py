import os
from pathlib import Path

import numpy as np
from astropy.io import fits
from astroquery.esa.xmm_newton import XMMNewton
from astroquery.utils.tap.core import TapPlus
from loguru import logger
from pysas.wrapper import Wrapper as sas


def run_xsa_query_async(query: str):
    """
    Use this if you expect more than 2000 returned elements
    """
    tap = TapPlus(url="https://nxsa.esac.esa.int/tap-server/tap")
    job = tap.launch_job_async(query=query)
    res = job.get_results()
    tap.remove_jobs([job.jobid])
    return res


def run_xsa_query(query: str):
    return XMMNewton.query_xsa_tap(query=query)


def download_data(observation_id: str, extension: str, filename: str) -> None:
    XMMNewton.download_data(observation_id, extension=extension, filename=filename)


def check_pps_dir(pps_dir: Path):
    """

    PURPOSE:
        Check if the input folder `pps_dir` is indeed a folder with XMM-SAS pipeline products (PPS)

    INPUT:
        pps_dir - str,
            Relative or absolute path to the pps_dir

    METHOD:
        Will check if the following files are available in pps_dir, in this order
        1. If pps_dir exists and it is directory
        1. CALINDEX file *OBXCALIND*
        1. Flaring background files: *FBKTSR*
        1. Calibrated event lists: *IEVLI*
        1. Attitude HK data: *ATTTSR*

        If any of these checks fails then we return False, i.e. it's not a PPS folder

        In case all files are found, then will return a dict with the absolute paths to each of these:
        'ccf_file', 'att_file', 'fbk_files', 'evl_files'
    """
    # check1
    if not pps_dir.exists():
        raise NotADirectoryError(f"Directory {pps_dir} not found!")

    pps_files = {}

    # check 2
    ccf = list(pps_dir.glob("*CALIND*"))
    if len(ccf) < 1:
        raise FileNotFoundError(f"Calibration index file *CALIND* file not found in {pps_dir}.")
    pps_files["ccf_file"] = ccf[0]

    # check 3
    fbk = list(pps_dir.glob("*FBKTSR*"))
    if len(fbk) < 1:
        raise FileNotFoundError(f"Flaring background *FBKTSR* files not found in {pps_dir}.")
    pps_files["fbk_files"] = fbk

    # check 4
    evl = list(pps_dir.glob("*IEVLI*"))
    if len(evl) < 1:
        raise FileNotFoundError(f"Calibrated event lists *IEVLI* files not found in {pps_dir}.")
    pps_files["evl_files"] = evl

    # check 5
    att = list(pps_dir.glob("*ATTTSR*"))
    if len(att) < 1:
        raise FileNotFoundError(f"Attitude *ATTTSR* file not found in {pps_dir}.")
    pps_files["att_file"] = att[0]

    return pps_files


def make_gti_pps(
    pps_files: dict[str, Path | list[Path]],
    instrument="all",
    out_dir: Path = None,
    max_expo=-1.0,
    verbose=True,
):
    """
    PURPOSE:
        Generate good-time-interval file using the XMM Pipeline Produced (PPS) flaring background and the PPS threshold

    INPUTS:
        pps_dir - str,
            the folder with the PPS files, checks will be performed to make sure it is indeed PPS folder
        instrument - str,
            One of m1, m2, pn or all, default all, case-insensitive. Which instrument to process.
        out_dir - str,
            the output folder where the .gti files will be saved. Default current working dir
        max_expo - float,
            Limit the GTI to have maximum exposure of `max_expo` ks. Set it to negative to skip this.
        plot_it - bool,
            whether to plot the flaring background curve with the GTI threshold
        save_plot - str,
            Name of the plot file, if needed. If None then no saving.

    METHOD:
        Will search in `pps_dir` for the PPS flaring background file *FBKTSR*
        The file header should have a keyword `FLCUTTHR` with the value of the background threshold
        Will run `tabgtigen` on the FLBKTSR table with expression `RATE <= FLTCUTTHR` and save it to .gti file
        Optionally, a plot will be produced and saved to a file (also optional)
    """
    if out_dir is None:
        out_dir = Path.cwd()
    # pps_files will contain the absolute path to all necessary files
    # rate_lim = {'m1': 0.0, 'm2': 0.0, 'pn': 0.0}
    inst_short = {"EMOS1": "m1", "EMOS2": "m2", "EPN": "pn"}

    gti_names = []
    for j in pps_files["fbk_files"]:
        if (instrument.upper() in j.stem.upper()) or (instrument.upper() == "ALL"):
            with fits.open(j) as hdu:
                # skip RGS files
                if "FLCUTTHR" not in hdu["RATE"].header:
                    logger.info(f"Cannot find threshold FLCUTTHR in {j}, cannot filter for GTI")
                    continue

                inst = hdu[0].header["INSTRUME"]
                if inst not in inst_short:
                    continue

                rate_lim = hdu["RATE"].header["FLCUTTHR"]
            # run gtigen
            gti_name = out_dir / f"{inst_short[inst]}_pps.gti"
            args = [
                f"table={j}",
                f'expression="RATE<={rate_lim}"',
                f"gtiset={gti_name.absolute().as_posix()}",
                "mingtisize=50.0",
            ]

            logFile = "DEFAULT" if verbose else os.devnull

            sas("tabgtigen", args, logFile).run()

            if max_expo > 0.0:
                # filter with max exposure
                if verbose:
                    logger.debug(f"Now creating GTI with limited exposure of {max_expo} ks")
                xgti_name = out_dir / f"{inst_short[inst]}_pps_{max_expo:.1f}ks.gti"
                max_expo_gti(gti_name, xgti_name, max_expo=max_expo)
            else:
                xgti_name = gti_name

            with fits.open(xgti_name, mode="update") as hdu:
                hdu["STDGTI"].header["METHOD"] = (
                    "pps",
                    "Method used to derive the rate threshold",
                )
                hdu["STDGTI"].header["RLIM"] = (rate_lim, "The PPS derived threshold")

            gti_names.append(xgti_name)

    return gti_names


def max_expo_gti(gti_infile: Path, gti_outfile: Path, max_expo=10.0) -> None:
    """
    PURPOSE:
        Using an input GTI file, will modify it in such a way as to have GTI duration of max_expo

    INPUTS:
        gti_infile - str,
            The name of the input GTI file, e.g. created using flaring background time series and a given rate threshold
        gti_outfile - str,
            The name of the output GTI file, filtering by this file will produce event list with max_expo duration
        max-expo - float,
            The required maximum exposure in ks
    """
    if not gti_infile.exists():
        raise FileNotFoundError(f"Input GTI file {gti_infile} not found")

    max_expo_sec = max_expo * 1000.0

    with fits.open(gti_infile) as hdu:
        nrec = len(hdu["STDGTI"].data)
        mask = np.zeros(nrec, dtype=bool)
        delta_time = hdu["STDGTI"].data["STOP"] - hdu["STDGTI"].data["START"]
        # first the easiest, find if there are GTI greater or equal to max_expo
        ix = np.where(delta_time >= max_expo_sec)[0]
        if len(ix) == 1:
            mask[ix] = 1
            hdu["STDGTI"].data["STOP"][ix] = hdu["STDGTI"].data["START"][ix] + max_expo_sec
        elif len(ix) > 1:
            imax = np.argmax(delta_time)
            mask[imax] = 1
            hdu["STDGTI"].data["STOP"][imax] = hdu["STDGTI"].data["START"][imax] + max_expo_sec
        else:
            # no single GTI is larger than max_expo and we have to accumulate them, starting with the largest
            # starting from largest GTI
            ixsort = delta_time.argsort()[-nrec:][::-1]
            sum_gti = 0.0
            for js in ixsort:  # TODO
                sum_gti += delta_time[js]
                if sum_gti >= max_expo_sec:
                    # last GTI will have to make it to max_expo_sec
                    mask[js] = 1
                    dd = sum_gti - max_expo_sec
                    hdu["STDGTI"].data["STOP"][js] = hdu["STDGTI"].data["START"][js] + dd
                    break
                mask[js] = 1
        hdu["STDGTI"].data = hdu["STDGTI"].data[mask]
        hdu.writeto(gti_outfile, overwrite=True)


def filter_events_gti(
    event_list: Path,
    gti_file: Path,
    pps_files: dict[str, Path | list[Path]],
    w_dir: Path,
    output_name=None,
    filter_expression=None,
    verbose=False,
):
    """
    PURPOSE:
        Filter event list with a GTI file

    INPUTS:
        event_list - str,
            FITS file with event list to filter with GTI
        gti_file - str,
            GTI file with periods of good time intervals
        max_expo - float,
            The maximum exposure time to extract, in ks.
            If negative, then it will not limit the exposure and will use the full GTI-filtered exposure.
        output_name - str,
            The name of the filtered output event list
        verbose - bool,
            If to echo some verbose information

        METHOD:
            Applying make_gti_pps() using the PPS flared background and the PPS derived threshold to create a GTI file
            then use this GTI file to filter the event list

    """
    # check if the input files exist
    if not event_list.exists():
        raise FileNotFoundError(f"Input event list {event_list} not found.")

    if not gti_file.exists():
        raise FileNotFoundError(f"GTI file {gti_file} not found.")

    # read some keywords for housekeeping
    hdr = fits.getheader(event_list, "EVENTS")
    inst = hdr["INSTRUME"]
    ontime0 = hdr["ONTIME"]

    # read some keywords for housekeeping from the GTI file, method and threshold
    hdr2 = fits.getheader(gti_file, "STDGTI")
    method = hdr2["METHOD"]

    # mapping the instrument names
    xinst = {"EMOS1": "m1", "EMOS2": "m2", "EPN": "pn"}

    # now filtering the event lists with the GTI
    if filter_expression is None:
        if "EMOS" in inst:
            expr = f"#XMMEA_EM && gti({gti_file},TIME) && (PI>150) && (PATTERN <= 12)"
        elif "EPN" in inst:
            expr = f"(FLAG == 0) && gti({gti_file},TIME) && (PI>150) && (PATTERN <= 4)"
    else:
        # user provided filtering
        expr = f"{filter_expression} && gti({gti_file},TIME)"

    output_name = w_dir / f"{xinst[inst]}_{method}_clean.evl" if output_name is None else w_dir / output_name

    args = [
        f"table={event_list}",
        "withfilteredset=Y",
        f"filteredset={output_name.absolute().as_posix()}",
        "destruct=Y",
        "keepfilteroutput=Y",
        f'expression="{expr}"',
    ]

    logFile = "DEFAULT" if verbose else os.devnull

    sas("evselect", args, logFile).run()

    if verbose:
        hdr1 = fits.getheader(output_name, "EVENTS")
        ontime1 = hdr1["ONTIME"]
        print(
            f"Input event list on-time {ontime0:.1f} s,\n"
            f"filtered list on-time {ontime1:.1f} s,\n"
            f"good time fraction {(100 * ontime1 / ontime0):.1f} %"
        )
    return output_name


def make_detxy_image(
    event_list: Path,
    w_dir: Path,
    pps_dir: Path,
    pps_files: dict[str, Path | list[Path]],
    output_name=None,
    low_energy=500,
    high_energy=15000,
    bin_size=80,
    radec_image=True,
    verbose=False,
):
    """
    PURPOSE:
        Generate image in DETX, DETY and add proper sky WCS from an event list

    INPUTS:
        event_list - str,
            FITS file with events list
        pps_dir - str,
            The folder with PPS products
        output_name - str,
            The name of the output image, if None then will be called {inst}_detxy_image.fits
        low_energy - int,
            The low PI energy for the image in eV, default 500 eV
        high_energy - int,
            The high PI energy for the image in eV, default 2000 eV
        bin_size - int,
            The bin_size in units of 0.05 arcsec, default 80 => 4" pixel
        radec_image - bool,
            If an image in sky coordinates (X,Y) should be also produced (default yes)
        verbose - bool,
            If to echo some verbose information

    METHOD:
        1. Standard run with evselect and create image with DETX, DETY
        2. Run ecoordconv to fix a reference pixel and RA, Dec to make a proper sky WCS
    """

    # check if the input files exist
    if not event_list.exists():
        raise FileNotFoundError(f"Input event list {event_list} not found.")

    # figure out the instrument from the event list header
    hdr = fits.getheader(event_list, "EVENTS")
    inst = hdr["INSTRUME"]
    # mapping the instrument names
    xinst = {"EMOS1": "m1", "EMOS2": "m2", "EPN": "pn"}

    logger.debug(f"*** {xinst[inst]}: generating image in DETX,DETY in band [{low_energy},{high_energy}] eV")

    if output_name is None:
        output_name = w_dir / f"{xinst[inst]}_{low_energy}_{high_energy}_detxy_image.fits"
    else:
        output_name = w_dir / output_name

    if "M1" in inst or "M2" in inst:
        expr = f"PI in [{low_energy}:{high_energy}] &&  (FLAG & 0x766ba000)==0 && PATTERN in [0:12]"
    else:
        expr = f"PI in [{low_energy}:{high_energy}] &&  FLAG==0 && PATTERN in [0:4]"

    args = [
        f"table={event_list}",
        "xcolumn=DETX",
        "ycolumn=DETY",
        f"ximagebinsize={bin_size}",
        f"yimagebinsize={bin_size}",
        f'expression="{expr}"',
        f"imageset={output_name.absolute().as_posix()}",
    ]

    logFile = "DEFAULT" if verbose else os.devnull
    sas("evselect", args, logFile).run()

    if verbose:
        print(f"\t DETXY image {output_name} created")

    # for reference will do an image in Sky coordinates
    if radec_image:
        if verbose:
            logger.info(f"*** {xinst[inst]}: generating image in RA,DEC in band [{low_energy},{high_energy}] eV")
        image_name_radec = output_name.absolute().as_posix().replace("detxy", "radec")
        args = [
            f"table={event_list}",
            "xcolumn=X",
            "ycolumn=Y",
            f"ximagebinsize={bin_size}",
            f"yimagebinsize={bin_size}",
            "squarepixels=yes",
            f'expression="{expr}"',
            f"imageset={image_name_radec}",
        ]

        evselect = sas("evselect", args)
        evselect.run()

        if verbose:
            print(f"\t RADEC image {image_name_radec} created")

    # now let's try to add WCS
    # will use `ecoordconv` with DETX,DETY=0,0 and get the IMAGE_X, IMAGE_Y pixel coordinates
    # and the corresponding RA, DEC
    # if verbose:
    #     print("Running ecoordconv")
    # args = [
    #     f"imageset={output_name.absolute().as_posix()}", "x=0", "y=0", "coordtype=det"
    # ]
    # ecoordconv = sas("ecoordconv", args)
    # ecoordconv.run()

    # for iline in status.stdout.decode().split("\n"):
    #     if "IM_X:" in iline:
    #         q = iline.split()
    #         xima = q[2]
    #         yima = q[3]
    #     if "DEC:" in iline:
    #         q = iline.split()
    #         ra = q[2]
    #         dec = q[3]

    # if verbose:
    #     print(f"Update the header of {output_name} with a new WCS")

    # with fits.open(output_name, mode="update") as hdu:
    #     header = hdu[0].header
    #     # Create a new WCS object.  The number of axes must be set
    #     # from the start
    #     header["CRVAL1"] = float(ra)
    #     header["CRVAL2"] = float(dec)
    #     header["CRPIX1"] = float(xima)
    #     header["CRPIX2"] = float(yima)
    #     cdelt1 = bin_size * header["REFYCDLT"]
    #     cdelt2 = -bin_size * header["REFXCDLT"]
    #     header["CDELT1"] = cdelt1
    #     header["CDELT2"] = cdelt2
    #     header["CTYPE1"] = "RA---TAN"
    #     header["CTYPE2"] = "DEC--TAN"

    #     # rotation
    #     crota2 = 90.0 - float(header["PA_PNT"])
    #     header["CROT2"] = crota2
    #     crota2_rad = math.radians(crota2)
    #     # add the CD matrix, just in case?
    #     cd1_1 = cdelt1 * math.cos(crota2_rad)
    #     cd1_2 = -1.0 * cdelt2 * math.sin(crota2_rad)
    #     cd2_1 = cdelt1 * math.sin(crota2_rad)
    #     cd2_2 = cdelt2 * math.cos(crota2_rad)
    #     header["CD1_1"] = cd1_1
    #     header["CD1_2"] = cd1_2
    #     header["CD2_1"] = cd2_1
    #     header["CD2_2"] = cd2_2
    #     header["COMMENT"] = "WCS added by IvanV"

    return output_name
