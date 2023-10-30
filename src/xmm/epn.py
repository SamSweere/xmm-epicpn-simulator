from typing import Tuple

import numpy as np
from astropy.io import fits

from src.xmm.xmm_ccf import get_epn_lincoord
from src.xmm.xmm_ccf import get_xmm_miscdata


def get_img_width_height(res_mult: int = 1) -> Tuple[int, int]:
    xrval, yrval = np.absolute(get_xyrval())

    p_delt = get_pixel_size(res_mult)

    max_x = round(float(np.max(xrval)), 3)
    max_y = round(float(np.max(yrval)), 3)

    size_x = np.ceil((max_x * 2 + 64 * p_delt * res_mult) / p_delt)
    size_y = np.ceil((max_y * 2 + 200 * p_delt * res_mult) / p_delt)

    return int(size_x), int(size_y)


def get_surface(res_mult: int = 1) -> float:
    """
    Returns:
        float: The surface of EPIC-pn in mm²
    """
    pixel_size = get_pixel_size(res_mult=res_mult)
    width, height = get_img_width_height(res_mult=res_mult)

    return (pixel_size ** 2) * width * height


def get_ccd_width_height(res_mult: int = 1) -> Tuple[int, int]:
    return 64 * res_mult, 200 * res_mult


def get_cc12_txy() -> Tuple[float, float]:
    epn_lincoord = get_epn_lincoord()
    with fits.open(name=epn_lincoord, mode="readonly") as file:
        header = file[1].header
        cc12_tx = header["CC12_TX"]
        cc12_ty = header["CC12_TY"]
    return cc12_tx, cc12_ty


def get_arc_mm_xy() -> Tuple[float, float]:
    epn_lincoord = get_epn_lincoord()
    with fits.open(name=epn_lincoord, mode="readonly") as file:
        header = file[1].header
        arc_mm_x = header["ARC_MM_Y"]
        arc_mm_y = header["ARC_MM_X"]
    return arc_mm_x, arc_mm_y


def get_xyrval() -> Tuple[np.ndarray, np.ndarray]:
    epn_lincoord = get_epn_lincoord()
    with fits.open(name=epn_lincoord, mode="readonly") as file:
        lincoord = file[1].data
        xrval = lincoord["X0"].astype(float)
        yrval = lincoord["Y0"].astype(float)

    return xrval, yrval


def get_pixel_size(res_mult: int = 1) -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data  # First entry is a PrimaryHDU, which is irrelevant for us
        epn = miscdata[miscdata["INSTRUMENT_ID"] == "EPN"]
        p_delt = epn[epn["PARM_ID"] == "MM_PER_PIXEL_X"]["PARM_VAL"].astype(float).item()  # Size of one pixel

    return round(p_delt / res_mult, 3)


def get_cdelt(res_mult: int = 1) -> float:
    # cdelt give the pixel sizes in degrees
    # cdelt from XMM_MISCDATA_0022.CCF PLATE_SCALE_X, the unit is in arsec, arsec to degree by deciding it by 3600
    xmm_miscdata = get_xmm_miscdata()
    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data
        epn = miscdata[miscdata["INSTRUMENT_ID"] == "EPN"]
        c_delt = epn[epn["PARM_ID"] == "PLATE_SCALE_X"]["PARM_VAL"].astype(float).item()

    c_delt = round((c_delt / 3600) / res_mult, 6)

    return c_delt


def get_focal_length() -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data  # First entry is a PrimaryHDU, which is irrelevant for us
        xrt3 = miscdata[miscdata["INSTRUMENT_ID"] == "XRT3"]  # XRT3 is the telescope, where EPN is located
        focallength = xrt3[xrt3["PARM_ID"] == "FOCAL_LENGTH"]["PARM_VAL"].astype(float).item()

    return focallength


def get_fov() -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data  # First entry is a PrimaryHDU, which is irrelevant for us
        xrt3 = miscdata[miscdata["INSTRUMENT_ID"] == "XRT3"]  # XRT3 is the telescope, where EPN is located
        fov = xrt3[xrt3["PARM_ID"] == "FOV_RADIUS"]["PARM_VAL"].astype(float).item() * 2  # Notice the 'RADIUS'

    return fov


def create_detector_mask(res_mult: int = 1) -> np.ndarray:
    width, height = get_img_width_height(res_mult)
    pixel_size = get_pixel_size(res_mult)
    xrval, yrval = get_xyrval()

    mask = np.zeros((width, height))

    small_gap = int(np.ceil((round(float(yrval[1] - yrval[0]), 3) - (64 * pixel_size * res_mult)) / pixel_size))
    large_gap = int(np.ceil((round(float(yrval[0] - yrval[3]), 3) - (64 * pixel_size * res_mult)) / pixel_size))
    vertical_gap = int(np.ceil((round(float(xrval[0] - xrval[9]), 3) - (200 * pixel_size * res_mult)) / pixel_size))

    drows = int(np.ceil((round(float(yrval[9] - yrval[0]), 3) / pixel_size)))

    # --- Upper row ---
    start = drows
    end = start + 64 * res_mult
    for i in range(6):
        mask[start:end + 1, :(200 * res_mult + 1)] = 1
        gap = small_gap if i != 2 else large_gap
        start = end + gap
        end = start + 64 * res_mult

    # --- Bottom row ---
    start = 0
    end = start + 64 * res_mult
    for i in range(6):
        mask[start:end + 1, (200 * res_mult + vertical_gap):] = 1
        gap = small_gap if i != 2 else large_gap
        start = end + gap
        end = start + 64 * res_mult

    # TODO
    # For whatever reason the images in the fits files are flipped -> We also have to flip the detector mask
    # mask = np.flipud(np.fliplr(mask))

    return mask


# TODO Add get_bad_pixels


def get_crpix(res_mult: int = 1) -> Tuple[float, float]:
    width, height = get_img_width_height(res_mult)
    cc12tx, cc12ty = get_cc12_txy()
    p_delt = get_pixel_size(res_mult)
    shift_x = cc12tx / p_delt
    shift_y = cc12ty / p_delt
    xrpix = round(((width + 1) / 2.0) + shift_x, 6)
    yrpix = round(((height + 1) / 2.0) - shift_y, 6)
    return xrpix, yrpix
