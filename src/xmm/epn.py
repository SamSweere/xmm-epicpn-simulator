from pathlib import Path
from typing import Tuple

import numpy as np
from astropy.io import fits

from src.xmm.xmm_ccf import get_xmm_miscdata


# This is the version if you want to calculate the image size through EPN_LINCOORD.CCF
# Unfortunately, I'm not entirely sure if the resulting image size is correct.
# def get_img_width_height(res_mult: int = 1) -> Tuple[int, int]:
#     epn_lincoord = get_epn_lincoord()
#     with fits.open(name=epn_lincoord, mode="readonly") as file:
#         lincoord = file[1].data
#         xrval = lincoord["Y0"].astype(float)
#         yrval = -lincoord["X0"].astype(float)
#
#     p_delt = get_pixel_size(res_mult)
#
#     dy = round(yrval[2] - yrval[5], 3)
#     drows = round(yrval[9] - yrval[0], 3)
#     width = np.ceil((dy + 64 * p_delt * res_mult + drows) / p_delt)
#     dx = round(xrval[5] - xrval[8], 3)
#     height = np.ceil((dx + 200 * p_delt * res_mult) / p_delt)
#
#     return int(width), int(height)

def get_img_width_height(res_mult: int = 1) -> Tuple[int, int]:
    return 403 * res_mult, 411 * res_mult


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


# This is the version if you want to calculate the detector mask through EPN_LINCOORD.CCF
# Unfortunately, I'm not entirely sure if the resulting mask is correct.
# def create_detector_mask(res_mult: int = 1) -> np.ndarray:
#     width, height = get_img_width_height(res_mult)
#     mask = np.zeros((width, height))
#
#     pixel_size = get_pixel_size(res_mult)
#
#     epn_lincoord = get_epn_lincoord()
#     with fits.open(name=epn_lincoord, mode="readonly") as file:
#         lincoord = file[1].data
#         xrval = lincoord["Y0"].astype(float)
#         yrval = -lincoord["X0"].astype(float)
#
#     small_gap = int(np.ceil((round(yrval[1] - yrval[0], 3) - (64 * pixel_size * res_mult)) / pixel_size))
#     large_gap = int(np.ceil((round(yrval[0] - yrval[3], 3) - (64 * pixel_size * res_mult)) / pixel_size))
#     vertical_gap = int(np.ceil((round(xrval[0] - xrval[9], 3) - (200 * pixel_size * res_mult)) / pixel_size))
#
#     drows = int(np.ceil((round(yrval[9] - yrval[0], 3) / pixel_size)))
#
#     # --- Upper row ---
#     start = drows
#     end = start + 64 * res_mult
#     for i in range(6):
#         mask[start:end + 1, :(200 * res_mult + 1)] = 1
#         gap = small_gap if i != 2 else large_gap
#         start = end + gap
#         end = start + 64 * res_mult
#
#     # Lower row
#     start = 0
#     end = start + 64 * res_mult
#     for i in range(6):
#         mask[start:end + 1, (200 * res_mult + vertical_gap):] = 1
#         gap = small_gap if i != 2 else large_gap
#         start = end + gap
#         end = start + 64 * res_mult
#
#     # For whatever reason the images in the fits files are flipped -> We also have to flip the detector mask
#     # mask = np.flipud(np.fliplr(mask))
#
#     return mask


def create_detector_mask(out_path: Path, res_mult: int = 1) -> np.ndarray:
    with fits.open('pn_expo_500_2000_detxy.ds') as hdu:
        # now make it a mask image
        mask = hdu[0].data > 0.0
        if res_mult > 1:
            mask = mask.repeat(res_mult, axis=0).repeat(res_mult, axis=1)

        hdu_out = hdu.copy()
        hdu_out[0].data = mask.astype(np.ubyte)
        hdu_out[0].scale('ubyte')

        hdu_out.writeto(out_path / f'pn_mask_500_2000_detxy_{res_mult}x.ds', overwrite=True)
    return mask


def get_crpix(res_mult: int = 1) -> Tuple[float, float]:
    width, height = get_img_width_height(res_mult)
    return round(width / 2.0, 6), round(height / 2.0, 6)
