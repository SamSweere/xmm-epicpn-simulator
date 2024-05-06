from pathlib import Path
from typing import Literal

import numpy as np
from loguru import logger
from lxml.etree import Element, ElementTree, SubElement

from src.xmm.utils import get_psf_file, get_vignet_file


def create_pn_xml(
    out_dir: Path,
    res_mult: int,
    xmm_filter: Literal["thin", "med", "thick"],
    sim_separate_ccds: bool,
    wait_time: float = 23.04e-6,  # Setting this to 0.0 eliminates out of time events
) -> list[Path]:
    from src.xmm.epn import get_focal_length, get_fov, get_pixel_size

    focallength = get_focal_length()
    fov = get_fov()
    p_delt = get_pixel_size(res_mult)

    # Change units from mm to m
    # See: http://www.sternwarte.uni-erlangen.de/~sixte/data/simulator_manual.pdf
    # in chap. "C: XML Instrument Configuration"
    focallength = round(focallength * 1e-3, 6)
    p_delt = round(p_delt * 1e-3, 6)

    if sim_separate_ccds:
        from src.xmm.epn import get_cc12_txy, get_ccd_width_height, get_xyrval

        max_x, max_y = get_ccd_width_height(res_mult=res_mult)
        xrval, yrval = get_xyrval()
        cc12tx, cc12ty = get_cc12_txy()
        xrval = (xrval - cc12tx) * 1e-3
        yrval = (yrval - cc12ty) * 1e-3
    else:
        from src.xmm.epn import get_cc12_txy, get_naxis12

        max_x, max_y = get_naxis12(res_mult=res_mult)
        xrval, yrval = get_cc12_txy()
        xrval = np.asarray([-xrval * 1e-3])
        yrval = np.asarray([-yrval * 1e-3])

    xrpix = round((max_x + 1) / 2.0, 6)
    yrpix = round((max_y + 1) / 2.0, 6)

    xml_paths: list[Path] = []
    loops = 12 if sim_separate_ccds else 1
    for i in range(loops):
        instrument = Element("instrument", telescop="XMM", instrume="EPN")

        telescope = SubElement(instrument, "telescope")
        # Based on the pixel fov and the biggest axis
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
        detector = SubElement(instrument, "detector", type="EPN")
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
        SubElement(detector, "rmf", filename=f"pn-{xmm_filter}-10.rmf")
        SubElement(detector, "arf", filename=f"pn-{xmm_filter}-10.arf")
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
            xml_path = out_dir / f"ccd{i + 1:02d}.xml"
            tree.write(xml_path, encoding="UTF-8", xml_declaration=True, pretty_print=True)
        else:
            xml_path = out_dir / "combined.xml"
            tree.write(xml_path, encoding="UTF-8", xml_declaration=True, pretty_print=True)

        xml_paths.append(xml_path)
    return xml_paths


def get_pn_xml(
    xml_dir: Path,
    res_mult: int,
    xmm_filter: Literal["thin", "med", "thick"],
    sim_separate_ccds: bool,
) -> list[Path]:
    instrument_path = xml_dir / "epn"
    root = instrument_path / xmm_filter / f"{res_mult}x"

    glob_pattern = "ccd*.xml" if sim_separate_ccds else "combined.xml"
    xml_paths: list[Path] = list(root.glob(glob_pattern))

    if sim_separate_ccds and len(xml_paths) != 12:
        logger.warning(
            f"'sim_separate_ccds' is set to 'True', but I could find only {len(xml_paths)} of the 12 CCDs."
            f"I will simulate only the CCDs given in these files. If that was intentional, then you can "
            f"ignore this warning. Otherwise abort the execution, create all XML files and re-run the"
            f"code."
        )

    return xml_paths


def create_mos_xml(
    out_dir: Path,
    emos_num: Literal[1, 2],
    res_mult: int,
    xmm_filter: Literal["thin", "med", "thick"],
    sim_separate_ccds: bool,
    wait_time: float = 23.04e-6,  # Setting this to 0.0 eliminates out of time events
) -> list[Path]:
    from src.xmm.emos import get_focal_length, get_fov, get_pixel_size

    if emos_num not in [1, 2]:
        raise ValueError(f"emos_num has to be either '1' or '2', but got {emos_num}!")

    focallength = get_focal_length(emos_num=emos_num)
    fov = get_fov(emos_num=emos_num)
    p_delt = get_pixel_size(emos_num=emos_num, res_mult=res_mult)

    # Change units from mm to m
    # See: http://www.sternwarte.uni-erlangen.de/~sixte/data/simulator_manual.pdf
    # in chap. "C: XML Instrument Configuration"
    focallength = round(focallength * 1e-3, 6)
    p_delt = round(p_delt * 1e-3, 6)

    if sim_separate_ccds:
        from src.xmm.emos import get_ccd_width_height, get_xyrval

        width, height = get_ccd_width_height(res_mult=res_mult)
        xrval, yrval = get_xyrval(emos_num=emos_num)
        xrval = xrval * 1e-3
        yrval = yrval * 1e-3
    else:
        from src.xmm.emos import get_img_width_height

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


def get_mos_xml(
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
        xml_paths.clear()

    return xml_paths
