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

from src.tools import xsa
from src.tools.files import decompress_targz
from src.xmm.ccf import get_emos_lincoord, get_telescope, get_xmm_miscdata
from src.xmm.tools import get_psf_file, get_vignet_file


def get_img_width_height(emos_num: Literal[1, 2], res_mult: int = 1) -> tuple[int, int]:
    xrval, yrval = np.absolute(get_xyrval(emos_num))

    p_delt = get_pixel_size(emos_num, res_mult)

    max_x = round(float(np.max(xrval)), 3)
    max_y = round(float(np.max(yrval)), 3)

    size_x = np.floor((max_x * 2 + 600 * p_delt * res_mult) / p_delt)
    size_y = np.floor((max_y * 2 + 600 * p_delt * res_mult) / p_delt)

    return int(size_x), int(size_y)


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


def get_cdelt(res_mult: int = 1) -> tuple[float, float]:
    # cdelt give the pixel sizes in degrees
    # The correct way would be to use get_plate_scale_xy()
    # BUT: When creating the dataset based on real observations
    # one has to give a binSize. This binSize is calculated as
    # follows: binSize = plate_scale / 0.05
    # Since the binSize has to be an integer, we can't use
    # the plate_scale given by the CCF (1.1008) and to
    # make our lifes easier for higher resolution images
    # we choose to set plate_scale to 1.0, which results
    # in binSize = 20
    c_delt = np.round((1.0 / 3600) / res_mult, 6)

    return c_delt, -c_delt


def get_naxis12(emos_num: Literal[1, 2], res_mult: int = 1) -> tuple[int, int]:
    if res_mult not in [1, 2, 4]:
        raise NotImplementedError
    if emos_num == 1:
        # These are the sizes one gets when creating the detector mask
        # images from real observation 0935190401 with binSize = 20
        # Since we rotate the image to represent the orientation
        # of EMOS1 relative to EPN we have to switch the axis
        if res_mult == 1:
            return 2006, 1317
        if res_mult == 2:
            return 4012, 2634
        if res_mult == 4:
            return 8024, 5267

    if emos_num == 2:
        # These are the sizes one gets when creating the detector mask
        # images from real observation 0935190401 with binSize = 20
        # Since we rotate the image to represent the orientation
        # of EMOS2 relative to EPN we don't have to switch the axis
        if res_mult == 1:
            return 1993, 2003
        if res_mult == 2:
            return 3986, 4006
        if res_mult == 4:
            return 7971, 8012


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


def create_mask(
    emin: float,
    emax: float,
    emos_num: Literal[1, 2],
    observation_id: str,
    out_dir: Path,
    mask_level: str,
    res_mults: list[int] = None,
) -> dict[int, Path]:
    if res_mults is None:
        res_mults = [1]
    inst = f"M{emos_num}"
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
        gtis = xsa.make_gti_pps(files, out_dir=tmp_dir, verbose=False)

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

        masks = {}
        for res_mult in res_mults:
            bin_size = 20 / res_mult

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
            expimgset = f"{inst}_expmap_{res_mult}x.fits"
            args = [
                f"imageset={imageset.resolve()}",
                f"attitudeset={odf_dir / 'atthk.dat'}",
                f"eventset={filtered.resolve()}",
                f"expimageset={expimgset}",
                f"pimin={int(emin * 1000)}",
                f"pimax={int(emax * 1000)}",
                "withdetcoords=true",
            ]
            sas("eexpmap", args, os.devnull).run()
            # Create emask
            emask = f"{inst}_emask_{res_mult}x.fits"
            sas("emask", [f"expimageset={expimgset}", f"detmaskset={emask}"], os.devnull).run()

            # Move to out_dir
            if mask_level == "expmap":
                mask_path = Path(out_dir) / "expmap" / expimgset
                with fits.open(Path.cwd() / expimgset, mode="update") as f:
                    f[0].data[f[0].data > 0] = 1

            if mask_level == "emask":
                mask_path = Path(out_dir) / "emask" / emask

            mask_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(mask_path.name, mask_path)

            masks[res_mult] = mask_path

    os.chdir(old_cwd)
    return masks


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

        if emos_num == 1:
            rotas = ["0.0", "90.0", "90.0", "90.0", "270.0", "270.0", "270.0"]
            yrval = -yrval
        else:
            rotas = ["0.0", "270.0", "270.0", "270.0", "90.0", "90.0", "90.0"]
            # We need to switch the axis for EMOS2 to be orthogonal to EMOS1
            # In SIXTE the x-axis points north and y-axis points west.
            # See https://xmmweb.esac.esa.int/docs/documents/CAL-MAN-0001.pdf
            # on page 5 for a visualisation.
            xrval, yrval = yrval, xrval

    else:
        width, height = get_img_width_height(emos_num=emos_num, res_mult=res_mult)
        xrval = yrval = 0.0
        xrval = np.asarray([xrval * 1e-3])
        yrval = np.asarray([yrval * 1e-3])

        rotas = ["0.0"]

    xrpix = round((width + 1) / 2.0, 6)
    yrpix = round((height + 1) / 2.0, 6)

    instrument = Element("instrument", telescop="XMM", instrume=f"EM{emos_num}")

    telescope = SubElement(instrument, "telescope")
    SubElement(telescope, "rmf", filename=f"mos{emos_num}-{xmm_filter}-10.rmf")
    SubElement(telescope, "arf", filename=f"mos{emos_num}-{xmm_filter}-10.arf")
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

    for i in range(len(rotas)):
        if (emos_num == 1) and (i == 2 or i == 5):
            # TODO Make the choice if CCD3 and CCD6 for EMOS1 should be used a config parameter
            continue
        detector = SubElement(instrument, "detector", type="ccd", chip=f"{i}")
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
        xml_path = out_dir / f"seperate_ccds_{xmm_filter}.xml"
    else:
        xml_path = out_dir / f"combined_ccd_{xmm_filter}.xml"

    tree.write(xml_path, encoding="UTF-8", xml_declaration=True, pretty_print=True)

    return xml_path


def get_xml(
    xml_dir: Path,
    emos_num: Literal[1, 2],
    res_mult: int,
    xmm_filter: Literal["thin", "med", "thick"],
    sim_separate_ccds: bool,
) -> Path:
    instrument_path = xml_dir / f"emos{emos_num}"
    root = instrument_path / xmm_filter / f"{res_mult}x"

    glob_pattern = f"seperate_ccds_{xmm_filter}.xml" if sim_separate_ccds else f"combined_ccd_{xmm_filter}.xml"
    xml_path: Path = next(root.glob(glob_pattern))

    if not xml_path:
        raise FileNotFoundError(f"Couldn't find {glob_pattern} for EMOS{emos_num} in {root.resolve()}!")

    return xml_path
