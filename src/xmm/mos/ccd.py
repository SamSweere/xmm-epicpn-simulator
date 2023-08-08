import json
from pathlib import Path
from lxml.etree import Element, SubElement, tostring, parse, ElementTree
from astropy.io import fits


class MosCcd:
    def __init__(
            self,
            ccf_path: Path,
            mos_num: int,
            num: int,
            res_mult: int,
            individual_ccd: bool,
            rotation
    ):
        if individual_ccd:
            raise NotImplementedError("Simulation of individual MOS-CCDs has not been implemented yet."
                                      "I can't find the size of an individual MOS-CCD anywhere")
        # TODO Add parameter which gives if we want a single CCD or the whole camera
        # TODO Add parameter to distinguish between EMOS1 and EMOS2
        # TODO If EMOS1, black out the dead ccds?
        if not float(res_mult).is_integer():
            error_message = "Resolution multiplication (res_mult) has to be a whole number"
            raise Exception(error_message)

        if not ccf_path.exists():
            raise NotADirectoryError(f"Could not find the ccf directory at '{ccf_path.resolve()}'!")

        if not 0 < mos_num < 3:
            raise ValueError(f"'mos_num' has to be either 1 or 2")

        if individual_ccd & (not 0 < num < 8):
            raise ValueError(f"'num' has to be an integer between 1 and 7!")

        mos_lincoord = list(ccf_path.glob(f"EMOS{mos_num}_LINCOORD*.CCF"))
        if not mos_lincoord:
            raise FileExistsError(f"Could not find any EPN_LINCOORD*.CCF in '{ccf_path.resolve()}'!")
        mos_lincoord.sort()
        mos_lincoord = mos_lincoord[-1]

        with fits.open(name=mos_lincoord, mode="readonly") as file:
            # See: https://xmmweb.esac.esa.int/docs/documents/CAL-MAN-0001.pdf chap. 4.3.20
            # CC12_TX/CC12_TY describe the transformation from CAMCOORD1 to CAMCOORD2 (in mm).
            # The origin in CAMCOORD1 is the geometrical center of the camera (center of chip wafer),
            # in CAMCOORD2 it is the point where the optical axis of the telescope in front the camera
            # intersects the focal plane.
            # Since xval and yval are locations in the focal plane, we need to shift them by cc12_tx and cc12_ty
            # respectively. (See: http://www.sternwarte.uni-erlangen.de/~sixte/data/simulator_manual.pdf chap.
            # "C: XML Instrument Configuration").
            lincoord = file[1].data
            if individual_ccd:
                xrval = lincoord["X0"][num - 1]
                yrval = lincoord["Y0"][num - 1]
            else:
                xrval = 0
                yrval = 0
        print(f"Finished reading {mos_lincoord.resolve()}")

        xmm_miscdata = list(ccf_path.glob("XMM_MISCDATA*.CCF"))
        if not xmm_miscdata:
            raise FileExistsError(f"Could not find any XMM_MISCDATA*.CCF in '{ccf_path.resolve()}'!")
        xmm_miscdata.sort()
        xmm_miscdata = xmm_miscdata[-1]

        with fits.open(name=xmm_miscdata, mode="readonly") as file:
            miscdata = file[1].data  # First entry is a PrimaryHDU, which is irrelevant for us

            # XRT1 = Telescope with EMOS1
            # XRT2 = Telescope with EMOS2
            xrt = miscdata[miscdata["INSTRUMENT_ID"] == f"XRT{mos_num}"]
            focallength = xrt[xrt["PARM_ID"] == "FOCAL_LENGTH"]["PARM_VAL"].astype(float).item()
            fov = xrt[xrt["PARM_ID"] == "FOV_RADIUS"]["PARM_VAL"].astype(float).item() * 2  # Notice the 'RADIUS'
            fov = fov + 0.3  # Add some margin. The 0.5 fov created too small images.

            emos = miscdata[miscdata["INSTRUMENT_ID"] == f"EMOS{mos_num}"]
            # Size of one pixel
            p_delt = emos[emos["PARM_ID"] == "MM_PER_PIXEL_X"]["PARM_VAL"].astype(float).item()
        print(f"Finished reading {xmm_miscdata.resolve()}")

        # TODO Also allow loading an .xml file instead of creation
        self.num = num

        self.instrument = Element("instrument", telescope="XMM", instrument="EPIC-MOS")

        self.telescope = SubElement(self.instrument, "telescope")
        # Based on the pixel fov and the biggest axis
        SubElement(self.telescope, "focallength", value=f"{focallength}")
        SubElement(self.telescope, "fov", diameter=f"{fov}")
        # TODO Check if the psf file exists
        # TODO Create psf files
        SubElement(self.telescope, "psf", filename=f"pn_psf_small_str_{1.0 / res_mult}_e_0.5_2.0_kev_v.1.3.fits")

        # Change units from mm to m
        # See: http://www.sternwarte.uni-erlangen.de/~sixte/data/simulator_manual.pdf
        # in chap. "C: XML Instrument Configuration"
        # Round the number in order to not have floating point errors
        # TODO Sam has added two pixels to xrval for the bottom ccds. Why?
        xrval = round(xrval * 1e-3 * res_mult, 12)
        yrval = round(yrval * 1e-3 * res_mult, 12)
        p_delt = round((p_delt * 1e-3) / res_mult, 12)
        # TODO Find source
        if individual_ccd:
            width = 64
            height = 200
        else:
            width = 600
            height = 600
        # The coordinates of the sensor are always in the middle of the sensor
        # Round the numbers in order to not have small floating point errors in the xml files
        xrpix = round((width * res_mult) / 2.0, 12) + 1
        yrpix = round((height * res_mult) / 2.0, 12) + 1

        # TODO I don't get the rota key. Why did Sam set this to 90 for ccds 1-6 and to 270 for ccds 7-12?
        rota = round(rotation, 12)
        self.detector = SubElement(self.instrument, 'detector', type='EPIC-pn')
        SubElement(self.detector, 'dimensions', xwidth=f"{width}", ywidth=f"{height}")
        SubElement(self.detector, 'wcs', xrpix=f"{xrpix}", yrpix=f"{yrpix}", xrval=f"{xrval}", yrval=f"{yrval}",
                   xdelt=f"{p_delt}", ydelt=f"{p_delt}", rota=f"{rota}")
        SubElement(self.detector, 'cte', value="1")
        # TODO Check if the rmf file exists
        SubElement(self.detector, 'rmf', filename="epn_e3_ff20_sdY9_v19.0.rmf")
        # TODO Check if the arf file exists
        SubElement(self.detector, 'arf', filename="pn_thin_onaxis.arf")
        # TODO: one could add instrument noise here, this is currently incorperated in the background files
        # if self.g_settings['noise']:
        #     # SubElement(detector, 'phabackground', filename = "pntffu_background_spectrum.fits")
        #     pass
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

    def pprint_xml(self):
        print(tostring(self.instrument, encoding="unicode", pretty_print=True))

    def save_xml(self, out_path: Path, prefix: str = ""):
        out_file = out_path / f"{prefix}_ccd{self.num:02d}.xml"
        self.tree.write(out_file, encoding='UTF-8', xml_declaration=True, pretty_print=True)
