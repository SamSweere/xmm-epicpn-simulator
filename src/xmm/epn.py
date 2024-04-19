import numpy as np
from astropy.io import fits

from src.xmm.ccf import get_epn_lincoord, get_telescope, get_xmm_miscdata


def get_max_xy(res_mult: int = 1) -> tuple[int, int]:
    xrval, yrval = np.absolute(get_xyrval())

    p_delt = get_pixel_size(res_mult)

    max_x = round(float(np.max(xrval)), 3)
    max_y = round(float(np.max(yrval)), 3)

    size_x = np.ceil((max_x * 2 + 64 * p_delt * res_mult) / p_delt)
    size_y = np.ceil((max_y * 2 + 200 * p_delt * res_mult) / p_delt)

    return int(size_x), int(size_y)


def get_naxis12(res_mult: int = 1) -> tuple[int, int]:
    fov_deg = get_fov()
    arc_mm_x, arc_mm_y = get_plate_scale_xy()

    fov_arcsec = fov_deg * 3600
    naxis1 = int(np.ceil(fov_arcsec / arc_mm_x))
    naxis2 = int(np.ceil(fov_arcsec / arc_mm_y))

    return naxis1 * res_mult, naxis2 * res_mult


def get_surface(res_mult: int = 1) -> float:
    """
    Returns:
        float: The surface of EPIC-pn in mm²
    """
    pixel_size = get_pixel_size(res_mult=res_mult)
    x, y = get_max_xy(res_mult=res_mult)

    return (pixel_size**2) * x * y


def get_ccd_width_height(res_mult: int = 1) -> tuple[int, int]:
    return 64 * res_mult, 200 * res_mult


def get_cc12_txy() -> tuple[float, float]:
    epn_lincoord = get_epn_lincoord()
    with fits.open(name=epn_lincoord, mode="readonly") as file:
        header = file[1].header
        cc12_tx = header["CC12_TX"]
        cc12_ty = header["CC12_TY"]
    return cc12_tx, cc12_ty


def get_plate_scale_xy() -> tuple[float, float]:
    xmm_miscdata = get_xmm_miscdata()
    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data
        epn = miscdata[miscdata["INSTRUMENT_ID"] == "EPN"]
        plate_scale_x = epn[epn["PARM_ID"] == "PLATE_SCALE_X"]["PARM_VAL"].astype(float).item()
        plate_scale_y = epn[epn["PARM_ID"] == "PLATE_SCALE_Y"]["PARM_VAL"].astype(float).item()
    return plate_scale_x, plate_scale_y


def get_xyrval() -> tuple[np.ndarray, np.ndarray]:
    epn_lincoord = get_epn_lincoord()
    with fits.open(name=epn_lincoord, mode="readonly") as file:
        lincoord = file[1].data
        xrval = lincoord["X0"].astype(float)
        yrval = lincoord["Y0"].astype(float)

    return xrval, yrval


def get_pixel_size(res_mult: int = 1) -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        # First entry is a PrimaryHDU, which is irrelevant for us
        miscdata = file[1].data
        epn = miscdata[miscdata["INSTRUMENT_ID"] == "EPN"]
        # Size of one pixel
        p_delt = epn[epn["PARM_ID"] == "MM_PER_PIXEL_X"]["PARM_VAL"].astype(float).item()

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
        # First entry is a PrimaryHDU, which is irrelevant for us
        miscdata = file[1].data
        telescope = get_telescope("epn")
        xrt = miscdata[miscdata["INSTRUMENT_ID"] == telescope]
        focallength = xrt[xrt["PARM_ID"] == "FOCAL_LENGTH"]["PARM_VAL"].astype(float).item()

    return focallength


def get_fov() -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        # First entry is a PrimaryHDU, which is irrelevant for us
        miscdata = file[1].data
        telescope = get_telescope("epn")
        xrt = miscdata[miscdata["INSTRUMENT_ID"] == telescope]
        # Notice the 'RADIUS'
        fov = xrt[xrt["PARM_ID"] == "FOV_RADIUS"]["PARM_VAL"].astype(float).item() * 2

    return fov


# WARNING: THE FOLLOWING FUNCTION WILL CREATE A THEORETICAL MASK FOR THE EPIC-PN DETECTOR
# FOR REALISM USE THE ACTUAL CALIBRATED DETECTOR MASK
def create_detector_mask(res_mult: int = 1) -> np.ndarray:
    width, height = get_max_xy(res_mult)
    pixel_size = get_pixel_size(res_mult)
    xrval, yrval = get_xyrval()

    mask = np.zeros((width, height))

    small_gap = int(np.ceil((round(float(yrval[1] - yrval[0]), 3) - (64 * pixel_size * res_mult)) / pixel_size))
    large_gap = int(np.ceil((round(float(yrval[0] - yrval[3]), 3) - (64 * pixel_size * res_mult)) / pixel_size))
    vertical_gap = int(np.ceil((round(float(xrval[0] - xrval[9]), 3) - (200 * pixel_size * res_mult)) / pixel_size))

    drows = int(np.ceil(round(float(yrval[9] - yrval[0]), 3) / pixel_size))

    # --- Upper row ---
    start = drows
    end = start + 64 * res_mult
    for i in range(6):
        mask[start : end + 1, : (200 * res_mult + 1)] = 1
        gap = small_gap if i != 2 else large_gap
        start = end + gap
        end = start + 64 * res_mult

    # --- Bottom row ---
    start = 0
    end = start + 64 * res_mult
    for i in range(6):
        mask[start : end + 1, (200 * res_mult + vertical_gap) :] = 1
        gap = small_gap if i != 2 else large_gap
        start = end + gap
        end = start + 64 * res_mult

    # TODO
    # For whatever reason the images in the fits files are flipped -> We also have to flip the detector mask
    # mask = np.flipud(np.fliplr(mask))

    return mask


# TODO Add get_bad_pixels


def get_shift_xy(res_mult: int = 1) -> tuple[float, float]:
    cc12tx, cc12ty = get_cc12_txy()
    p_delt = get_pixel_size(res_mult)
    shift_x = cc12tx / p_delt
    shift_y = cc12ty / p_delt
    return shift_x, shift_y
