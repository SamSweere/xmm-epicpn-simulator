from pathlib import Path
from typing import Optional

from lxml.etree import Element, ElementTree, SubElement, tostring

from utils.log import plog, elog


class XmmPnCcd:
    def __init__(self, num, g_settings, xrval, yrval, rotation, verbose):
        self.num = num
        self.g_settings = g_settings
        self.verbose = verbose
        self.name = g_settings['base_name'] + "_ccd{:02d}".format(num)
        self.file_location = ""
        self.tree: Optional[ElementTree] = None

        # Based on the pixel fov and the biggest axis
        self.fov_diameter = g_settings['fov']
        self.focal_length = g_settings['focal_length']

        self.res_mult = g_settings['res_mult']
        if not float(self.res_mult).is_integer():
            error_message = "Resolution multiplication (res_mult) has to be a whole number"
            elog(error_message)
            raise Exception(error_message)

        self.width = g_settings['width']
        self.height = g_settings['height']
        # The coorinates of the sensor are always in the middle of the sensor
        # Round the numbers in order to not have small floating point errors in the xml files
        self.xrpix = round(self.width / 2.0, 12) + 1
        self.yrpix = round(self.height / 2.0, 12) + 1
        self.xrval = round(xrval,
                           12)  # Round the number in order to not have floating point errors
        self.yrval = round(yrval, 12)
        self.p_delt = round(g_settings['p_delt'], 12)
        self.rota = round(rotation, 12)

    def _create_tree(self) -> None:
        instrument = Element('instrument', telescop="XMM", instrument="EPIC-pn")
        telescope = SubElement(instrument, 'telescope')
        SubElement(telescope, 'focallength', value=str(self.focal_length * 1.0))
        SubElement(telescope, 'fov', diameter=str(self.fov_diameter * 1.0))
        if self.res_mult == 1:
            # Use the psf for he normal resolution
            # SubElement(telescope, 'psf', filename='pn_psf_str_1.0_v.1.3.fits')
            SubElement(telescope, 'psf', filename='pn_psf_small_str_1.0_e_0.5_2.0_kev_v.1.3.fits')
        elif self.res_mult == 2:
            # SubElement(telescope, 'psf', filename='pn_psf_str_0.5_v.1.3.fits')
            SubElement(telescope, 'psf', filename='pn_psf_small_str_0.5_e_0.5_2.0_kev_v.1.3.fits')
        elif self.res_mult == 4:
            # SubElement(telescope, 'psf', filename='pn_psf_str_0.25_v.1.3.fits')
            SubElement(telescope, 'psf', filename='pn_psf_small_str_0.25_e_0.5_2.0_kev_v.1.3.fits')
        else:
            raise ValueError(f"There is no psf file for the res_mult: {self.res_mult}")

        detector = SubElement(instrument, 'detector', type='EPIC-pn')

        SubElement(detector, 'dimensions', xwidth=str(self.width), ywidth=str(self.height))
        SubElement(detector, 'wcs', xrpix=str(self.xrpix * 1.0), yrpix=str(self.yrpix * 1.0),
                   xrval=str(self.xrval * 1.0),
                   yrval=str(self.yrval * 1.0), xdelt=str(self.p_delt * 1.0), ydelt=str(self.p_delt * 1.0),
                   rota=str(self.rota * 1.0))
        SubElement(detector, 'cte', value="1")

        SubElement(detector, 'rmf', filename="epn_e3_ff20_sdY9_v19.0.rmf")

        SubElement(detector, 'arf', filename="pn_thin_onaxis.arf")

        # TODO: one could add instrument noise here, this is currently incorperated in the background files
        # if self.g_settings['noise']:
        #     # SubElement(detector, 'phabackground', filename = "pntffu_background_spectrum.fits")
        #     pass

        SubElement(detector, 'vignetting', filename="xmm_pn_vignet_0_6_fov.fits")
        SubElement(detector, 'split', type="gauss", par1=str(11.e-6 / self.res_mult))
        SubElement(detector, 'threshold_readout_lo_keV', value="0.")
        SubElement(detector, 'threshold_event_lo_keV', value="200.e-3")
        SubElement(detector, 'threshold_split_lo_fraction', value="0.01")
        SubElement(detector, 'threshold_pattern_up_keV', value="12.")

        readout = SubElement(detector, 'readout', mode="time")
        SubElement(readout, 'wait', time=str(self.g_settings['readout_wait_time']))

        loop = SubElement(readout, 'loop', start="0", end=str(self.height - 1), increment="1", variable="$i")
        SubElement(loop, 'readoutline', lineindex="0", readoutindex="$i")
        SubElement(loop, 'lineshift')
        SubElement(loop, 'wait', time=str(self.g_settings['lineshift_wait_time']))

        SubElement(readout, 'newframe')

        self.tree = ElementTree(instrument)

    def get_tree(self) -> ElementTree:
        # TODO Check why the return type hint is not propagated
        if self.tree is None:
            self._create_tree()

        return self.tree

    def pprint_xml(self):
        print(tostring(self.get_tree(), encoding="unicode", pretty_print=True))

    def save_xml(self, folder_location: Path):
        self.file_location = folder_location / f"{self.name}.xml"
        self.get_tree().write(self.file_location, encoding='UTF-8', xml_declaration=True, pretty_print=True)
        plog("Saved ccd{:02d}".format(self.num) + " to: " + self.file_location, verbose=self.verbose)
