from pathlib import Path

from astropy.io import fits
from lxml.etree import Element, SubElement, tostring, parse, ElementTree


class PnCcd:
    def __init__(
            self,
            ccf_path: Path,
            num: int,
            res_mult: int,
            individual_ccd: bool,
            rotation: int,
            noise: bool = False
    ):
        # TODO Add parameter which gives if we want a single CCD or the whole camera
        if not float(res_mult).is_integer():
            error_message = "Resolution multiplication (res_mult) has to be a whole number"
            raise Exception(error_message)

        if not ccf_path.exists():
            raise NotADirectoryError(f"Could not find the ccf directory at '{ccf_path.resolve()}'!")

        if individual_ccd & (not 0 < num < 13):
            raise ValueError(f"'num' has to be an integer between 1 and 12!")

        epn_lincoord = list(ccf_path.glob("EPN_LINCOORD*.CCF"))
        if not epn_lincoord:
            raise FileExistsError(f"Could not find any EPN_LINCOORD*.CCF in '{ccf_path.resolve()}'!")
        epn_lincoord.sort()
        epn_lincoord = epn_lincoord[-1]

        with fits.open(name=epn_lincoord, mode="readonly") as file:
            # See: https://xmmweb.esac.esa.int/docs/documents/CAL-MAN-0001.pdf chap. 4.3.20
            # CC12_TX/CC12_TY describe the transformation from CAMCOORD1 to CAMCOORD2 (in mm).
            # The origin in CAMCOORD1 is the geometrical center of the camera (center of chip wafer),
            # in CAMCOORD2 it is the point where the optical axis of the telescope in front the camera
            # intersects the focal plane.
            # Since xval and yval are locations in the focal plane, we need to shift them by cc12_tx and cc12_ty
            # respectively. (See: http://www.sternwarte.uni-erlangen.de/~sixte/data/simulator_manual.pdf chap.
            # "C: XML Instrument Configuration").
            cc12_tx = file[1].header["CC12_TX"]
            cc12_ty = file[1].header["CC12_TY"]
            lincoord = file[1].data
            if individual_ccd:
                # TODO Why did Sam rotate these values?
                xrval = lincoord["X0"][num - 1] + cc12_tx
                yrval = lincoord["Y0"][num - 1] + cc12_ty
            else:
                xrval = cc12_tx
                yrval = cc12_ty
        print(f"Finished reading {epn_lincoord.resolve()}")

        xmm_miscdata = list(ccf_path.glob("XMM_MISCDATA*.CCF"))
        if not xmm_miscdata:
            raise FileExistsError(f"Could not find any XMM_MISCDATA*.CCF in '{ccf_path.resolve()}'!")
        xmm_miscdata.sort()
        xmm_miscdata = xmm_miscdata[-1]

        with fits.open(name=xmm_miscdata, mode="readonly") as file:
            miscdata = file[1].data  # First entry is a PrimaryHDU, which is irrelevant for us

            xrt3 = miscdata[miscdata["INSTRUMENT_ID"] == "XRT3"]  # XRT3 is the telescope, where EPN is located
            focallength = xrt3[xrt3["PARM_ID"] == "FOCAL_LENGTH"]["PARM_VAL"].astype(float).item()
            fov = xrt3[xrt3["PARM_ID"] == "FOV_RADIUS"]["PARM_VAL"].astype(float).item() * 2  # Notice the 'RADIUS'
            fov = fov + 0.3  # Add some margin. The 0.5 fov created too small images.

            epn = miscdata[miscdata["INSTRUMENT_ID"] == "EPN"]
            p_delt = epn[epn["PARM_ID"] == "MM_PER_PIXEL_X"]["PARM_VAL"].astype(float).item()  # Size of one pixel
        print(f"Finished reading {xmm_miscdata.resolve()}")

        # TODO Also allow loading an .xml file instead of creation
        self.num = num

        # Change units from mm to m
        # See: http://www.sternwarte.uni-erlangen.de/~sixte/data/simulator_manual.pdf
        # in chap. "C: XML Instrument Configuration"
        # Round the number in order to not have floating point errors
        focallength = round(focallength * 1e-3, 12)
        # TODO Sam has added two pixels to xrval for the bottom ccds. Why?
        xrval = round(xrval * 1e-3 * res_mult, 12)
        yrval = round(yrval * 1e-3 * res_mult, 12)
        p_delt = round((p_delt * 1e-3) / res_mult, 12)
        # TODO Find source
        if individual_ccd:
            width = 64 * res_mult
            height = 200 * res_mult
        else:
            width = 403 * res_mult
            height = 411 * res_mult
        # The coordinates of the sensor are always in the middle of the sensor
        # Round the numbers in order to not have small floating point errors in the xml files
        xrpix = round(width / 2.0, 12) + 1
        yrpix = round(height / 2.0, 12) + 1

        # TODO I don't get the rota key. Why did Sam set this to 90 for ccds 1-6 and to 270 for ccds 7-12?
        rota = round(rotation, 12)

        self.instrument = Element("instrument", telescope="XMM", instrument="EPIC-MOS")

        self.telescope = SubElement(self.instrument, "telescope")
        # Based on the pixel fov and the biggest axis
        SubElement(self.telescope, "focallength", value=f"{focallength}")
        SubElement(self.telescope, "fov", diameter=f"{fov}")
        # TODO Check if the psf file exists
        SubElement(self.telescope, "psf", filename=f"m1_psf_small_str_{1.0 / res_mult}_e_0.5_2.0_kev_v.1.3.fits")
        self.detector = SubElement(self.instrument, 'detector', type='EPIC-MOS')
        SubElement(self.detector, 'dimensions', xwidth=f"{width}", ywidth=f"{height}")
        SubElement(self.detector, 'wcs', xrpix=f"{xrpix}", yrpix=f"{yrpix}", xrval=f"{xrval}", yrval=f"{yrval}",
                   xdelt=f"{p_delt}", ydelt=f"{p_delt}", rota=f"{rota}")
        SubElement(self.detector, 'cte', value="1")
        # TODO Check if the rmf file exists
        SubElement(self.detector, 'rmf', filename="m1_e19_im_pall_o.rmf")
        # TODO Check if the arf file exists
        SubElement(self.detector, 'arf', filename="mos1-thin-10.arf")
        if noise:
            SubElement(self.detector, 'phabackground', filename="pntffu_background_spectrum.fits")
        # TODO Check if the vignetting file exists
        SubElement(self.detector, 'vignetting', filename="xmm_pn_vignet_0_6_fov.fits")
        SubElement(self.detector, 'split', type="gauss", par1=f"{11.e-6 / res_mult}")
        SubElement(self.detector, 'threshold_readout_lo_keV', value="0.")
        SubElement(self.detector, 'threshold_event_lo_keV', value="200.e-3")
        SubElement(self.detector, 'threshold_split_lo_fraction', value="0.01")
        SubElement(self.detector, 'threshold_pattern_up_keV', value="12.")

        self.readout = SubElement(self.detector, 'readout', mode="time")
        SubElement(self.readout, 'wait', time="68.75e-3")

        self.loop = SubElement(self.readout, 'loop', start="0", end=f"{height - 1}", increment="1", variable="$i")
        SubElement(self.loop, 'readoutline', lineindex="0", readoutindex="$i")
        SubElement(self.loop, 'lineshift')
        SubElement(self.loop, 'wait', time="0.0")  # Setting this to 0.0 eliminates out of time events

        SubElement(self.readout, 'newframe')

        self.tree = ElementTree(self.instrument)

    @classmethod
    def from_file(cls, num: int, file: Path):
        tree = parse(file)
        instrument = tree.getroot()

        pass

    def get_focallength(self) -> float:
        return float(self.telescope[0].get("value"))

    def get_fov(self) -> float:
        return float(self.telescope[1].get("fov"))

    def get_psf(self) -> str:  # TODO Change to return the path to psf file
        return self.telescope[2].get("filename")

    def pprint_xml(self):
        print(tostring(self.instrument, encoding="unicode", pretty_print=True))

    def save_xml(self, out_path: Path, prefix: str = ""):
        out_file = out_path / f"{prefix}_ccd{self.num:02d}.xml"
        self.tree.write(out_file, encoding='UTF-8', xml_declaration=True, pretty_print=True)
