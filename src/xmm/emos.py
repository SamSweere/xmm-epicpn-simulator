from pathlib import Path
from typing import Literal

import numpy as np
from astropy.io import fits
from loguru import logger
from lxml.etree import Element, ElementTree, SubElement

from src.xmm.ccf import get_emos_lincoord, get_telescope, get_xmm_miscdata
from src.xmm.utils import get_psf_file, get_vignet_file


def get_img_width_height(emos_num: Literal[1, 2], res_mult: int = 1) -> tuple[int, int]:
    xrval, yrval = np.absolute(get_xyrval(emos_num))

    p_delt = get_pixel_size(emos_num, res_mult)

    max_x = round(float(np.max(xrval)), 3)
    max_y = round(float(np.max(yrval)), 3)

    size_x_mm = max_x * 2 + 600 * p_delt
    size_y_mm = max_y * 2 + 600 * p_delt
    
    size_x_px = np.ceil(size_x_mm / p_delt) * res_mult
    size_y_px = np.ceil(size_y_mm / p_delt) * res_mult

    return int(size_x_px), int(size_y_px)


def get_surface(emos_num: Literal[1, 2], res_mult: int = 1) -> float:
    pixel_size = get_pixel_size(emos_num, res_mult)
    width, height = get_img_width_height(emos_num, res_mult)

    return (pixel_size**2) * width * height


def get_ccd_width_height(res_mult: int = 1) -> tuple[int, int]:
    """
    Returns:
        Tuple[int, int]: CCD width and height in pixels.
    """
    return 600 * res_mult, 600 * res_mult


def get_plate_scale_xy(emos_num: Literal[1, 2]) -> tuple[float, float]:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data
        emos = miscdata[miscdata["INSTRUMENT_ID"] == f"EMOS{emos_num}"]
        plate_scale_x = emos[emos["PARM_ID"] == "PLATE_SCALE_X"]["PARM_VAL"].astype(float).item()
        plate_scale_y = emos[emos["PARM_ID"] == "PLATE_SCALE_Y"]["PARM_VAL"].astype(float).item()

    return plate_scale_x, plate_scale_y


def get_xyrval(emos_num: Literal[1, 2]) -> tuple[np.ndarray, np.ndarray]:
    emos_lincoord = get_emos_lincoord(emos_num=emos_num)
    with fits.open(name=emos_lincoord, mode="readonly") as file:
        lincoord = file[1].data
        # Use the primary readout node and not the redundant one.
        lincoord = lincoord[lincoord["NODE_ID"] == 0]
        xrval = lincoord["X0"].astype(float)
        yrval = lincoord["Y0"].astype(float)

    return xrval, yrval


def get_pixel_size(emos_num: Literal[1, 2], res_mult: int = 1) -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data
        emos = miscdata[miscdata["INSTRUMENT_ID"] == f"EMOS{emos_num}"]
        # Size of one pixel
        p_delt = emos[emos["PARM_ID"] == "MM_PER_PIXEL_X"]["PARM_VAL"].astype(float).item()

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
        # First entry is a PrimaryHDU, which is irrelevant for us
        miscdata = file[1].data
        telescope = get_telescope(f"emos{emos_num}")
        xrt = miscdata[miscdata["INSTRUMENT_ID"] == telescope]
        focallength = xrt[xrt["PARM_ID"] == "FOCAL_LENGTH"]["PARM_VAL"].astype(float).item()

    return focallength


def get_fov(emos_num: Literal[1, 2]) -> float:
    xmm_miscdata = get_xmm_miscdata()

    with fits.open(name=xmm_miscdata, mode="readonly") as file:
        miscdata = file[1].data
        telescope = get_telescope(f"emos{emos_num}")
        xrt = miscdata[miscdata["INSTRUMENT_ID"] == telescope]
        # Notice the 'RADIUS'
        fov = xrt[xrt["PARM_ID"] == "FOV_RADIUS"]["PARM_VAL"].astype(float).item() * 2

    return fov


def create_xml(
    out_dir: Path,
    emos_num: Literal[1, 2],
    res_mult: int,
    xmm_filter: Literal["thin", "med", "thick"],
    sim_separate_ccds: bool,
    wait_time: float = 23.04e-6,  # Setting this to 0.0 eliminates out of time events
) -> list[Path]:
    # Change units from mm to m
    # See: http://www.sternwarte.uni-erlangen.de/~sixte/data/simulator_manual.pdf
    # in chap. "C: XML Instrument Configuration"
    focallength = round(get_focal_length(emos_num=emos_num) * 1e-3, 6)
    p_delt = round(get_pixel_size(emos_num=emos_num, res_mult=res_mult) * 1e-3, 6)
    fov = get_fov(emos_num=emos_num)

    if sim_separate_ccds:
        width, height = get_ccd_width_height(res_mult=res_mult)
        xrval, yrval = get_xyrval(emos_num=emos_num)
        xrval = xrval * 1e-3
        yrval = yrval * 1e-3
    else:
        width, height = get_img_width_height(emos_num=emos_num, res_mult=res_mult)
        xrval = yrval = 0.0
        xrval = np.asarray([xrval * 1e-3])
        yrval = np.asarray([yrval * 1e-3])

    xrpix = round((width + 1) / 2.0, 6)
    yrpix = round((height + 1) / 2.0, 6)

    xml_paths: list[Path] = []
    if sim_separate_ccds:
        loops = 7
        if emos_num == 1:
            rotas = ["0.0", "90.0", "90.0", "90.0", "270.0", "270.0", "270.0"]
        else:
            rotas = ["270.0", "0.0", "0.0", "0.0", "180.0", "180.0", "180.0"]
    else:
        loops = 1
        rotas = ["0.0"]

    for i in range(loops):
        instrument = Element("instrument", telescop="XMM", instrume=f"EM{emos_num}")

        telescope = SubElement(instrument, "telescope")
        # Based on the pixel fov and the biggest axis
        SubElement(telescope, "focallength", value=f"{focallength}")
        SubElement(telescope, "fov", diameter=f"{fov}")
        SubElement(
            telescope,
            "psf",
            filename=f"{get_psf_file(xml_dir=out_dir, instrument_name=f'emos{emos_num}', res_mult=res_mult).name}",
        )
        SubElement(
            telescope,
            "vignetting",
            filename=f"{get_vignet_file(xml_dir=out_dir, instrument_name=f'emos{emos_num}').name}",
        )
        detector = SubElement(instrument, "detector", type=f"EM{emos_num}")
        SubElement(detector, "dimensions", xwidth=f"{width}", ywidth=f"{height}")
        SubElement(
            detector,
            "wcs",
            xrpix=f"{xrpix}",
            yrpix=f"{yrpix}",
            xrval=np.format_float_positional(xrval[i], 6),
            yrval=np.format_float_positional(yrval[i], 6),
            xdelt=f"{p_delt}",
            ydelt=f"{p_delt}",
            rota=f"{rotas[i]}",
        )
        SubElement(detector, "cte", value="1")
        SubElement(detector, "rmf", filename=f"mos{emos_num}-{xmm_filter}-10.rmf")
        SubElement(detector, "arf", filename=f"mos{emos_num}-{xmm_filter}-10.arf")
        SubElement(detector, "split", type="gauss", par1=f"{11.e-6 / res_mult}")
        SubElement(detector, "threshold_readout_lo_keV", value="0.")
        SubElement(detector, "threshold_event_lo_keV", value="200.e-3")
        SubElement(detector, "threshold_split_lo_fraction", value="0.01")
        SubElement(detector, "threshold_pattern_up_keV", value="15.")

        readout = SubElement(detector, "readout", mode="time")
        SubElement(readout, "wait", time="2.6")

        loop = SubElement(
            readout,
            "loop",
            start="0",
            end=f"{height - 1}",
            increment="1",
            variable="$i",
        )
        SubElement(loop, "readoutline", lineindex="0", readoutindex="$i")
        SubElement(loop, "lineshift")
        if sim_separate_ccds:
            SubElement(loop, "wait", time=f"{wait_time}")  # Setting this to 0.0 eliminates out of time events

        SubElement(readout, "newframe")

        tree = ElementTree(instrument)
        if sim_separate_ccds:
            xml_path = out_dir / f"ccd{i + 1:02d}.xml"
            tree.write(xml_path, encoding="UTF-8", xml_declaration=True, pretty_print=True)
        else:
            xml_path = out_dir / "combined.xml"
            tree.write(xml_path, encoding="UTF-8", xml_declaration=True, pretty_print=True)

        xml_paths.append(xml_path)
    return xml_paths


def get_xml(
    xml_dir: Path,
    emos_num: Literal[1, 2],
    res_mult: int,
    xmm_filter: Literal["thin", "med", "thick"],
    sim_separate_ccds: bool,
) -> list[Path]:
    instrument_path = xml_dir / f"emos{emos_num}"
    root = instrument_path / xmm_filter / f"{res_mult}x"

    glob_pattern = "ccd*.xml" if sim_separate_ccds else "combined.xml"
    xml_paths: list[Path] = list(root.glob(glob_pattern))

    if sim_separate_ccds and len(xml_paths) != 7:
        logger.warning(
            f"'sim_separate_ccds' is set to 'True', but I could find only {len(xml_paths)} of the 7 CCDs."
            f"I will simulate only the CCDs given in these files. If that was intentional, then you can "
            f"ignore this warning. Otherwise abort the execution, create all XML files and re-run the"
            f"code."
        )

    return xml_paths
