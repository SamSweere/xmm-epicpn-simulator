from typing import Literal, Tuple

import numpy as np
from astropy.io import fits

from src.xmm.ccf import get_xmm_miscdata, get_emos_lincoord


def get_img_width_height(emos_num: Literal[1, 2], res_mult: int = 1) -> Tuple[int, int]:
    xrval, yrval = np.absolute(get_xyrval(emos_num))

    p_delt = get_pixel_size(emos_num, res_mult)

    max_x = round(float(np.max(xrval)), 3)
    max_y = round(float(np.max(yrval)), 3)

    size_x = np.ceil((max_x * 2 + 600 * p_delt * res_mult) / p_delt)
    size_y = np.ceil((max_y * 2 + 600 * p_delt * res_mult) / p_delt)

    if emos_num == 1:
        return int(size_y), int(size_x)
    else:
        return int(size_x), int(size_y)


def get_surface(emos_num: Literal[1, 2], res_mult: int = 1) -> float:
    pixel_size = get_pixel_size(emos_num, res_mult)
    width, height = get_img_width_height(emos_num, res_mult)

    return (pixel_size ** 2) * width * height


def get_ccd_width_height(res_mult: int = 1) -> Tuple[int, int]:
    """
    Returns:
        Tuple[int, int]: CCD width and height in pixels.
    """
    return 600 * res_mult, 600 * res_mult


def get_cc12_txy() -> Tuple[float, float]:
    return 0.0, 0.0


def get_xyrval(emos_num: Literal[1, 2]) -> Tuple[np.ndarray, np.ndarray]:
    emos_lincoord = get_emos_lincoord(emos_num=emos_num)
    with fits.open(name=emos_lincoord, mode="readonly") as file:
        lincoord = file[1].data
        lincoord = lincoord[lincoord["NODE_ID"] == 0]  # Use the primary readout node and not the redundant one.
        xrval = lincoord["X0"].astype(float)
        yrval = lincoord["Y0"].astype(float)

    return xrval, yrval


def get_pixel_size(emos_num: Literal[1, 2], res_mult: int = 1) -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data
        emos = miscdata[miscdata["INSTRUMENT_ID"] == f"EMOS{emos_num}"]
        p_delt = emos[emos["PARM_ID"] == "MM_PER_PIXEL_X"]["PARM_VAL"].astype(float).item()  # Size of one pixel

    return round(p_delt / res_mult, 3)


def get_cdelt(emos_num: Literal[1, 2], res_mult: int = 1) -> float:
    # cdelt give the pixel sizes in degrees
    # cdelt from XMM_MISCDATA_0022.CCF PLATE_SCALE_X, the unit is in arsec, arsec to degree by deciding it by 3600
    xmm_miscdata = get_xmm_miscdata()
    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data
        emos = miscdata[miscdata["INSTRUMENT_ID"] == f"EMOS{emos_num}"]
        c_delt = emos[emos["PARM_ID"] == "PLATE_SCALE_X"]["PARM_VAL"].astype(float).item()

    c_delt = round((c_delt / 3600) / res_mult, 6)

    return c_delt


def get_focal_length(emos_num: Literal[1, 2]) -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data  # First entry is a PrimaryHDU, which is irrelevant for us
        xrt = miscdata[miscdata["INSTRUMENT_ID"] == f"XRT{emos_num}"]  # XRT3 is the telescope, where EPN is located
        focallength = xrt[xrt["PARM_ID"] == "FOCAL_LENGTH"]["PARM_VAL"].astype(float).item()

    return focallength


def get_fov(emos_num: Literal[1, 2]) -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data
        xrt = miscdata[miscdata["INSTRUMENT_ID"] == f"XRT{emos_num}"]
        fov = xrt[xrt["PARM_ID"] == "FOV_RADIUS"]["PARM_VAL"].astype(float).item() * 2  # Notice the 'RADIUS'

    return fov


def get_crpix(emos_num: Literal[1, 2], res_mult: int = 1) -> Tuple[float, float]:
    width, height = get_img_width_height(emos_num, res_mult)
    xrval, yrval = get_xyrval(emos_num)
    p_delt = get_pixel_size(emos_num, res_mult)
    shift_x = float(xrval[0]) / p_delt
    shift_y = float(yrval[0]) / p_delt
    xrpix = round(((width + 1) / 2.0) + shift_x, 6)
    yrpix = round(((height + 1) / 2.0) - shift_y, 6)
    return xrpix, yrpix
