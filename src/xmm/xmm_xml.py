import os
from pathlib import Path
from typing import Literal, List

import numpy as np
from astropy.io import fits
from lxml.etree import Element, SubElement, ElementTree

from src.xmm.epn import get_img_width_height, get_focal_length, get_fov, get_pixel_size
from src.xmm.xmm_ccf import get_epn_lincoord


def create_pn_xml(
        res_mult: int,
        xmm_filter: Literal["thin", "med", "thick"],
        sim_separate_ccds: bool,
        wait_time: float = 23.04e-6  # Setting this to 0.0 eliminates out of time events
) -> List[Path]:
    instrument_path = Path(os.environ["SIXTE"]) / "share" / "sixte" / "instruments" / "xmm" / "epicpn"
    if not instrument_path.exists():
        raise NotADirectoryError(f"It looks like you haven't downloaded the instrument files provided by SIXTE "
                                 f"(see https://www.sternwarte.uni-erlangen.de/sixte/instruments/)! Please download "
                                 f"and extract them as given in their instructions.")

    epn_lincoord = get_epn_lincoord()

    with fits.open(name=epn_lincoord, mode="readonly") as file:
        # See: https://xmmweb.esac.esa.int/docs/documents/CAL-MAN-0001.pdf chap. 4.3.20
        # See: http://www.sternwarte.uni-erlangen.de/~sixte/data/simulator_manual.pdf chap.
        # "C: XML Instrument Configuration".
        lincoord = file[1].data
        # I don't know why, but the axis have to be switched and the signs flipped for yrval
        xrval = lincoord["Y0"].astype(float)
        yrval = -lincoord["X0"].astype(float)
    print(f"Finished reading {epn_lincoord.resolve()}")

    focallength = get_focal_length()
    fov = get_fov()
    p_delt = get_pixel_size(res_mult)

    # Change units from mm to m
    # See: http://www.sternwarte.uni-erlangen.de/~sixte/data/simulator_manual.pdf
    # in chap. "C: XML Instrument Configuration"
    focallength = round(focallength * 1e-3, 6)
    xrval = xrval * 1e-3
    yrval = yrval * 1e-3
    p_delt = round((p_delt * 1e-3) / res_mult, 6)

    if sim_separate_ccds:
        width = 64 * res_mult
        height = 200 * res_mult
    else:
        width, height = get_img_width_height(res_mult)

    # The coordinates of the sensor are always in the middle of the sensor
    xrpix = round(width / 2.0, 6)
    yrpix = round(height / 2.0, 6)

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
        SubElement(detector, 'wcs', xrpix=f"{xrpix}", yrpix=f"{yrpix}",
                   xrval=np.format_float_positional(xrval[i], 6),
                   yrval=np.format_float_positional(yrval[i], 6),
                   xdelt=f"{p_delt}", ydelt=f"{p_delt}", rota=f"90.0")
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
            xml_path = instrument_path / f"ccd{i + 1:02d}_{xmm_filter}_{res_mult}x.xml"
            tree.write(xml_path, encoding='UTF-8', xml_declaration=True, pretty_print=True)
        else:
            xml_path = instrument_path / f"combined_{xmm_filter}_{res_mult}x.xml"
            tree.write(xml_path, encoding='UTF-8', xml_declaration=True, pretty_print=True)

        xml_paths.append(xml_path)
    return xml_paths


def get_pn_xml(
        res_mult: int,
        xmm_filter: Literal["thin", "med", "thick"],
        sim_separate_ccds: bool,
) -> List[Path]:
    instrument_path = Path(os.environ["SIXTE"]) / "share" / "sixte" / "instruments" / "xmm" / "epicpn"

    xml_paths: List[Path] = []
    if sim_separate_ccds:
        for i in range(12):
            xml_path = instrument_path / f"ccd{i + 1:02d}_{xmm_filter}_{res_mult}x.xml"
            if xml_path.exists():
                xml_paths.append(xml_path)
        if not len(xml_paths) == 12:
            xml_paths.clear()
    else:
        xml_path = instrument_path / f"combined_{xmm_filter}_{res_mult}x.xml"
        if xml_path.exists():
            xml_paths.append(xml_path)

    return xml_paths


def create_mos_xml(

):
    pass


if __name__ == '__main__':
    create_pn_xml(1, "thin", True)
