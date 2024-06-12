from pathlib import Path
from typing import Literal

import numpy as np
from astropy.io import fits
from lxml.etree import Element, ElementTree, SubElement

from src.xmm.ccf import get_epn_lincoord, get_telescope, get_xmm_miscdata
from src.xmm.utils import get_psf_file, get_vignet_file


def get_img_width_height(res_mult: int = 1) -> tuple[int, int]:
    xrval, yrval = np.absolute(get_xyrval())

    p_delt = get_pixel_size(res_mult)

    max_x = np.round(np.max(xrval), 3)
    max_y = np.round(np.max(yrval), 3)

    size_x = np.floor((max_x * 2 + 64 * p_delt * res_mult) / p_delt)
    size_y = np.floor((max_y * 2 + 200 * p_delt * res_mult) / p_delt)

    return int(size_x), int(size_y)


def get_surface(res_mult: int = 1) -> float:
    """
    Returns:
        float: The surface of EPIC-pn in mm²
    """
    pixel_size = get_pixel_size(res_mult=res_mult)
    x, y = get_img_width_height(res_mult=res_mult)

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
        plate_scale_x = epn[epn["PARM_ID"] == "PLATE_SCALE_X"]["PARM_VAL"].item()
        plate_scale_y = epn[epn["PARM_ID"] == "PLATE_SCALE_Y"]["PARM_VAL"].item()
    return plate_scale_x, plate_scale_y


def get_xyrval() -> tuple[np.ndarray, np.ndarray]:
    epn_lincoord = get_epn_lincoord()
    with fits.open(name=epn_lincoord, mode="readonly") as file:
        lincoord = file[1].data
        # TODO This is a temporary fix.
        # The two rows should be perfectly aligned on the x-axis,
        # but for whatever reason they are not. XMM-Newton Helpdesk
        # has been contacted. Answer is still pending.
        # This is the "correct" version:
        xrval = lincoord["X0"]
        yrval = lincoord["Y0"]
        # The following steps can be deleted when the issue is fixed
        xrval[-6:] = -xrval[-6:]

    return xrval, yrval


def get_pixel_size(res_mult: int = 1) -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        # First entry is a PrimaryHDU, which is irrelevant for us
        miscdata = file[1].data
        epn = miscdata[miscdata["INSTRUMENT_ID"] == "EPN"]
        # Size of one pixel
        p_delt = epn[epn["PARM_ID"] == "MM_PER_PIXEL_X"]["PARM_VAL"].item()

    return round(p_delt / res_mult, 3)


def get_cdelt(res_mult: int = 1) -> float:
    # cdelt give the pixel sizes in degrees
    # The correct way would be to use get_plate_scale_xy()
    # BUT: When creating the dataset based on real observations
    # one has to give a binSize. This binSize is calculated as
    # follows: binSize = plate_scale / 0.05
    # Since the binSize has to be an integer, we can't use
    # the plate_scale given by the CCF (4.12838) and to
    # make our lifes easier for higher resolution images
    # we choose to set plate_scale to 4.0, which results
    # in binSize = 80
    c_delt = np.round((4.0 / 3600) / res_mult, 6)

    return c_delt


def get_naxis12(res_mult: int = 1) -> tuple[int, int]:
    # 403x411 is the image size one gets when creating
    # images from real observations with binSize = 80
    return (403 * res_mult, 411 * res_mult)


def get_focal_length() -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        # First entry is a PrimaryHDU, which is irrelevant for us
        miscdata = file[1].data
        telescope = get_telescope("epn")
        xrt = miscdata[miscdata["INSTRUMENT_ID"] == telescope]
        focallength = xrt[xrt["PARM_ID"] == "FOCAL_LENGTH"]["PARM_VAL"].item()

    return focallength


def get_fov() -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        # First entry is a PrimaryHDU, which is irrelevant for us
        miscdata = file[1].data
        telescope = get_telescope("epn")
        xrt = miscdata[miscdata["INSTRUMENT_ID"] == telescope]
        # Notice the 'RADIUS'
        fov = xrt[xrt["PARM_ID"] == "FOV_RADIUS"]["PARM_VAL"].item() * 2

    return fov


# WARNING: THE FOLLOWING FUNCTION WILL CREATE A THEORETICAL MASK FOR THE EPIC-PN DETECTOR
# FOR REALISM USE THE ACTUAL CALIBRATED DETECTOR MASK
def create_detector_mask(res_mult: int = 1) -> np.ndarray:
    width, height = get_img_width_height(res_mult)
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


def create_xml(
    out_dir: Path,
    res_mult: int,
    xmm_filter: Literal["thin", "med", "thick"],
    sim_separate_ccds: bool,
    wait_time: float = 23.04e-6,  # Setting this to 0.0 eliminates out of time events
) -> list[Path]:
    # Change units from mm to m
    # See: http://www.sternwarte.uni-erlangen.de/~sixte/data/simulator_manual.pdf
    # in chap. "C: XML Instrument Configuration"
    focallength = round(get_focal_length() * 1e-3, 6)
    p_delt = round(get_pixel_size(res_mult) * 1e-3, 6)
    fov = get_fov()

    if sim_separate_ccds:
        max_y, max_x = get_ccd_width_height(res_mult=res_mult)
        yrval, xrval = get_xyrval()
        cc12tx, cc12ty = get_cc12_txy()
        xrval = np.round(xrval + cc12tx, 3) * 1e-3
        yrval = np.round(yrval - cc12ty, 3) * 1e-3
    else:
        max_x, max_y = get_img_width_height(res_mult=res_mult)
        xrval, yrval = get_cc12_txy()
        xrval = np.asarray([-xrval * 1e-3])
        yrval = np.asarray([-yrval * 1e-3])

    xrpix = round((max_x + 1) / 2.0, 6)
    yrpix = round((max_y + 1) / 2.0, 6)

    instrument = Element("instrument", telescop="XMM", instrume="EPN")

    telescope = SubElement(instrument, "telescope")
    SubElement(telescope, "rmf", filename=f"pn-{xmm_filter}-10.rmf")
    SubElement(telescope, "arf", filename=f"pn-{xmm_filter}-10.arf")
    SubElement(telescope, "focallength", value=f"{focallength}")
    SubElement(telescope, "fov", diameter=f"{fov}")
    SubElement(
        telescope,
        "psf",
        filename=f"{get_psf_file(xml_dir=out_dir, instrument_name='epn', res_mult=res_mult).name}",
    )
    SubElement(
        telescope,
        "vignetting",
        filename=f"{get_vignet_file(xml_dir=out_dir, instrument_name='epn').name}",
    )

    for i in range(len(xrval)):
        detector = SubElement(instrument, "detector", type="ccd", chip=f"{i}")
        SubElement(detector, "dimensions", xwidth=f"{max_x}", ywidth=f"{max_y}")
        # See https://www.aanda.org/articles/aa/pdf/2019/10/aa35978-19.pdf Appendix A about the rota
        SubElement(
            detector,
            "wcs",
            xrpix=f"{xrpix}",
            yrpix=f"{yrpix}",
            xrval=np.format_float_positional(xrval[i], 6),
            yrval=np.format_float_positional(yrval[i], 6),
            xdelt=f"{p_delt}",
            ydelt=f"{p_delt}",
            rota=f"{'180.0' if i < 6 else '0.0'}",
        )
        SubElement(detector, "cte", value="1")
        SubElement(detector, "split", type="gauss", par1=f"{11.e-6 / res_mult}")
        SubElement(detector, "threshold_readout_lo_keV", value="0.")
        SubElement(detector, "threshold_event_lo_keV", value="200.e-3")
        SubElement(detector, "threshold_split_lo_fraction", value="0.01")
        SubElement(detector, "threshold_pattern_up_keV", value="15.")

        readout = SubElement(detector, "readout", mode="time")
        SubElement(readout, "wait", time="68.75e-3")

        loop = SubElement(readout, "loop", start="0", end=f"{max_y - 1}", increment="1", variable="$i")
        SubElement(loop, "readoutline", lineindex="0", readoutindex="$i")
        SubElement(loop, "lineshift")
        if sim_separate_ccds:
            SubElement(loop, "wait", time=f"{wait_time}")  # Setting this to 0.0 eliminates out of time events

        SubElement(readout, "newframe")

    tree = ElementTree(instrument)

    if sim_separate_ccds:
        xml_path = out_dir / f"seperate_ccds_{xmm_filter}.xml"
    else:
        xml_path = out_dir / f"combined_ccd_{xmm_filter}.xml"

    tree.write(xml_path, encoding="UTF-8", xml_declaration=True, pretty_print=True)

    return xml_path


def get_xml(
    xml_dir: Path,
    res_mult: int,
    xmm_filter: Literal["thin", "med", "thick"],
    sim_separate_ccds: bool,
) -> Path:
    instrument_path = xml_dir / "epn"
    root = instrument_path / xmm_filter / f"{res_mult}x"

    glob_pattern = f"seperate_ccds_{xmm_filter}.xml" if sim_separate_ccds else "combined.xml"
    xml_path: Path = next(root.glob(glob_pattern))

    return xml_path
