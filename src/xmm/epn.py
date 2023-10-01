from typing import Tuple

import numpy as np
from astropy.io import fits

from src.xmm.xmm_ccf import get_epn_lincoord, get_xmm_miscdata


def get_img_width_height(res_mult: int = 1) -> Tuple[float, float]:
    epn_lincoord = get_epn_lincoord()
    with fits.open(name=epn_lincoord, mode="readonly") as file:
        lincoord = file[1].data
        xrval = lincoord["Y0"].astype(float)
        yrval = -lincoord["X0"].astype(float)

    p_delt = get_pixel_size(res_mult)

    dy = round(yrval[2] - yrval[5], 3)
    drows = round(yrval[9] - yrval[0], 3)
    width = np.ceil((dy + 64 * p_delt * res_mult + drows) / p_delt)
    dx = round(xrval[5] - xrval[8], 3)
    height = np.ceil((dx + 200 * p_delt * res_mult) / p_delt)

    return width, height


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
