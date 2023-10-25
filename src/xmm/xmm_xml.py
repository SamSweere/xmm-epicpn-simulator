import os
from pathlib import Path
from typing import Literal, List

import numpy as np
from lxml.etree import Element, SubElement, ElementTree
from loguru import logger


def create_pn_xml(
        res_mult: int,
        xmm_filter: Literal["thin", "med", "thick"],
        sim_separate_ccds: bool,
        wait_time: float = 23.04e-6  # Setting this to 0.0 eliminates out of time events
) -> List[Path]:
    from src.xmm.epn import get_focal_length, get_fov, get_pixel_size

    instrument_path = Path(os.environ["SIXTE"]) / "share" / "sixte" / "instruments" / "xmm" / "epicpn"
    if not instrument_path.exists():
        raise NotADirectoryError(f"It looks like you haven't downloaded the instrument files provided by SIXTE "
                                 f"(see https://www.sternwarte.uni-erlangen.de/sixte/instruments/)! Please download "
                                 f"and extract them as given in their instructions.")

    out_dir = instrument_path / xmm_filter / f"{res_mult}x"
    out_dir.mkdir(parents=True, exist_ok=True)
    # Add symbolic links to used files
    files = [f"pn_psf_{1.0 / res_mult}x_e_0.5_2.0_kev.fits", "xmm_pn_vignet.fits", f"pn-{xmm_filter}-10.rmf",
             f"pn-{xmm_filter}-10.arf"]
    for file in files:
        tmp_link = out_dir / file
        if not tmp_link.exists():
            tmp_link.symlink_to(instrument_path / file)

    focallength = get_focal_length()
    fov = get_fov()
    p_delt = get_pixel_size(res_mult)

    # Change units from mm to m
    # See: http://www.sternwarte.uni-erlangen.de/~sixte/data/simulator_manual.pdf
    # in chap. "C: XML Instrument Configuration"
    focallength = round(focallength * 1e-3, 6)
    p_delt = round(p_delt * 1e-3, 6)

    if sim_separate_ccds:
        from src.xmm.epn import get_ccd_width_height, get_xyrval, get_cc12_txy
        width, height = get_ccd_width_height(res_mult=res_mult)
        xrval, yrval = get_xyrval()
        cc12tx, cc12ty = get_cc12_txy()
        xrval = (xrval - cc12tx) * 1e-3
        yrval = (yrval - cc12ty) * 1e-3
    else:
        from src.xmm.epn import get_img_width_height, get_cc12_txy
        width, height = get_img_width_height(res_mult=res_mult)
        xrval, yrval = get_cc12_txy()
        xrval = np.asarray([xrval * 1e-3])
        yrval = np.asarray([yrval * 1e-3])

    xrpix = round((width + 1) / 2.0, 6)
    yrpix = round((height + 1) / 2.0, 6)

    xml_paths: List[Path] = []
    loops = 12 if sim_separate_ccds else 1
    for i in range(loops):
        instrument = Element("instrument", telescop="XMM", instrume="EPIC-PN")

        telescope = SubElement(instrument, "telescope")
        # Based on the pixel fov and the biggest axis
        SubElement(telescope, "focallength", value=f"{focallength}")
        SubElement(telescope, "fov", diameter=f"{fov}")
        SubElement(telescope, "psf", filename=f"pn_psf_{1.0 / res_mult}x_e_0.5_2.0_kev.fits")
        SubElement(telescope, 'vignetting', filename="xmm_pn_vignet.fits")
        detector = SubElement(instrument, 'detector', type='EPIC-PN')
        SubElement(detector, 'dimensions', xwidth=f"{width}", ywidth=f"{height}")
        # See https://www.aanda.org/articles/aa/pdf/2019/10/aa35978-19.pdf Appendix A about the rota
        SubElement(detector, 'wcs', xrpix=f"{xrpix}", yrpix=f"{yrpix}",
                   xrval=np.format_float_positional(xrval[i], 6),
                   yrval=np.format_float_positional(yrval[i], 6),
                   xdelt=f"{p_delt}", ydelt=f"{p_delt}", rota=f"{'180.0' if i < 6 else '0.0'}")
        SubElement(detector, 'cte', value="1")
        SubElement(detector, 'rmf', filename=f"pn-{xmm_filter}-10.rmf")
        SubElement(detector, 'arf', filename=f"pn-{xmm_filter}-10.arf")
        SubElement(detector, 'split', type="gauss", par1=f"{11.e-6 / res_mult}")
        SubElement(detector, 'threshold_readout_lo_keV', value="0.")
        SubElement(detector, 'threshold_event_lo_keV', value="200.e-3")
        SubElement(detector, 'threshold_split_lo_fraction', value="0.01")
        SubElement(detector, 'threshold_pattern_up_keV', value="12.")

        readout = SubElement(detector, 'readout', mode="time")
        SubElement(readout, 'wait', time="68.75e-3")

        loop = SubElement(readout, 'loop', start="0", end=f"{height - 1}", increment="1", variable="$i")
        SubElement(loop, 'readoutline', lineindex="0", readoutindex="$i")
        SubElement(loop, 'lineshift')
        if sim_separate_ccds:
            SubElement(loop, 'wait', time=f"{wait_time}")  # Setting this to 0.0 eliminates out of time events

        SubElement(readout, 'newframe')

        tree = ElementTree(instrument)
        if sim_separate_ccds:
            xml_path = out_dir / f"ccd{i + 1:02d}.xml"
            tree.write(xml_path, encoding='UTF-8', xml_declaration=True, pretty_print=True)
        else:
            xml_path = out_dir / f"combined.xml"
            tree.write(xml_path, encoding='UTF-8', xml_declaration=True, pretty_print=True)

        xml_paths.append(xml_path)
    return xml_paths


def get_pn_xml(
        res_mult: int,
        xmm_filter: Literal["thin", "med", "thick"],
        sim_separate_ccds: bool
) -> List[Path]:
    instrument_path = Path(os.environ["SIXTE"]) / "share" / "sixte" / "instruments" / "xmm" / "epicpn"
    root = instrument_path / xmm_filter / f"{res_mult}x"

    glob_pattern = "ccd*.xml" if sim_separate_ccds else "combined.xml"
    xml_paths: List[Path] = list(root.glob(glob_pattern))

    if sim_separate_ccds and not len(xml_paths) == 12:
        logger.warning(f"'sim_separate_ccds' is set to 'True', but I could find only {len(xml_paths)} of the 12 CCDs."
                       f"I will simulate only the CCDs given in these files. If that was intentional, then you can "
                       f"ignore this warning. Otherwise abort the execution, create all XML files and re-run the"
                       f"code.")

    return xml_paths


def create_mos_xml(

):
    pass


if __name__ == '__main__':
    create_pn_xml(1, "thin", True)
