import os
import shutil
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal

import numpy as np
from astropy.io import fits
from lxml.etree import Element, ElementTree, SubElement
from pysas.wrapper import Wrapper as sas

from src import xsa
from src.xmm.ccf import get_epn_lincoord, get_telescope, get_xmm_miscdata
from src.xmm.utils import get_psf_file, get_vignet_file
from src.xmm_utils.file_utils import decompress_targz


def get_img_width_height(res_mult: int = 1) -> tuple[int, int]:
    xrval, yrval = np.absolute(get_xyrval())

    p_delt = get_pixel_size(res_mult)

    max_x = np.round(np.max(xrval), 3)
    max_y = np.round(np.max(yrval), 3)

    size_x = np.floor((max_x * 2 + 64 * p_delt * res_mult) / p_delt)
    size_y = np.floor((max_y * 2 + 200 * p_delt * res_mult) / p_delt)

    return int(size_x), int(size_y)


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
    # These are the sizes one gets when creating the detector mask
    # images from real observation 0935190401 with binSize = 80
    if res_mult == 1:
        return 403, 411
    if res_mult == 2:
        return 805, 822
    if res_mult == 4:
        return 1609, 1643

    raise NotImplementedError


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


def create_emask(observation_id: str, out_dir: Path, res_mults: list[int] = None) -> dict[int, Path]:
    if res_mults is None:
        res_mults = [1]
    inst = "PN"
    old_cwd = os.getcwd()
    with TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        obs_dir = tmp_dir / observation_id
        os.chdir(tmp_dir)
        xsa.download_data(observation_id, "FTZ", observation_id)
        with tarfile.open(f"{observation_id}.tar") as tar:
            tar.extractall()
        pps_dir = obs_dir / "pps"
        files = xsa.check_pps_dir(pps_dir)
        gtis = xsa.make_gti_pps(files, verbose=False)

        assert len(gtis) > 0

        evl = None
        for file in files["evl_files"]:
            if inst in file.stem.upper():
                evl = file
                break
        assert evl is not None

        gti = None
        for file in gtis:
            if inst in file.stem.upper():
                gti = file
                break
        assert gti is not None

        filtered = xsa.filter_events_gti(evl, gti, files, output_name=f"{observation_id}_cleaned.fits", w_dir=obs_dir)

        assert filtered

        # Create atthkset
        odf_dir = obs_dir / "odf"
        os.chdir(odf_dir)
        decompress_targz(odf_dir / f"{observation_id}.tar.gz", odf_dir)
        sas("odfingest", [f"odfdir={odf_dir}", "withodfdir=true"]).run()

        sum_file = next(odf_dir.glob("*SUM.SAS"))
        sas("atthkgen", ["-o", f"{sum_file}"]).run()

        emasks = {}
        for res_mult in res_mults:
            bin_size = 80 / res_mult

            imageset = xsa.make_detxy_image(
                filtered,
                pps_dir=pps_dir,
                pps_files=files,
                bin_size=bin_size,
                output_name=f"{observation_id}_detxy.fits",
                w_dir=obs_dir,
                radec_image=False,
            )

            os.chdir(obs_dir)
            # Create eexpmap
            expimgset = f"pn_expmap_{res_mult}x.fits"
            args = [
                f"imageset={imageset.resolve()}",
                f"attitudeset={odf_dir / 'atthk.dat'}",
                f"eventset={filtered.resolve()}",
                f"expimageset={expimgset}",
                "withdetcoords=true",
            ]
            sas("eexpmap", args, os.devnull).run()
            # Create emask
            emask = f"pn_emask_{res_mult}x.fits"
            sas("emask", [f"expimageset={expimgset}", f"detmaskset={emask}"], os.devnull).run()

            # Move to out_dir
            emask_path = Path(out_dir / "emask" / emask)
            emask_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(emask, emask_path)
            emasks[res_mult] = emask_path

    os.chdir(old_cwd)
    return emasks


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
        xrval = np.round(xrval + cc12ty, 3) * 1e-3
        yrval = np.round(yrval - cc12tx, 3) * 1e-3
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

    if not xml_path:
        raise FileNotFoundError(f"Couldn't find {glob_pattern} for EPN in {root.resolve()}!")

    return xml_path
