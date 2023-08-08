import json
from pathlib import Path
from typing import Dict

from src.xmm.pn.ccd import PnCcd


class XmmPnFf:
    def __init__(self, res_mult=1, noise=True, sim_separate_ccds=False, verbose=True):
        self.verbose = verbose
        self.base_name = "fullframe_thinfilter"

        # TODO Find reference:
        DEFAULT_IMG_WIDTH = 403
        DEFAULT_IMG_HEIGHT = 411

        # Independent of sim_separate_ccds:
        # cdelt give the pixel sizes in degrees
        # cdelt from XMM_MISCDATA_0022.CCF PLATE_SCALE_X, the unit is in arsec, arsec to degree by deciding it by 3600
        # The cdelt is scaled in order to have the same layout as the real images. I did not find the reason why this is needed
        # The scaling factor is determined by the pixel difference on the x-axis
        self.cdelt_x = 4.0 / 3600 / float(res_mult)
        self.cdelt_y = 4.0 / 3600 / float(res_mult)

        # The total width and height of the image
        self.total_pixel_width = DEFAULT_IMG_WIDTH * res_mult  # 420 #*res_mult #402 # np.ceil((self.xyrval['06'][1] - self.xyrval['03'][1])/p_delt + 2*0.5*ccd_pixel_width + 1) # add one for the two half pixels we are missing
        self.total_pixel_height = DEFAULT_IMG_HEIGHT * res_mult  # 420 #*res_mult #400 # np.ceil((self.xyrval['01'][0] - self.xyrval['10'][0])/p_delt + 2*0.5*ccd_pixel_height + 1)

        # CRPIX1 and CRPIX2 are the point corresponding to the optical axis on the instrument define the corresponding
        # pixel coordinate

        # Method based on knowing the boresight pixel on the final image
        # If the resolution of the image changes this will become incorrect
        crpix1_corr = 0
        crpix2_corr = 0

        if res_mult == 2:
            crpix1_corr += -0.5  # Seems to be correct
            crpix2_corr += -0.5  # Seems to be correct
        elif res_mult == 4:
            crpix1_corr += -1.5
            crpix2_corr += -1.5

        self.crpix1 = 244 * res_mult + crpix1_corr
        self.crpix2 = 224 * res_mult + crpix2_corr

        self.ccds = []
        # TODO Fix path
        ccf_path = Path("/home/bojantodorkov/Projects/xmm-epicpn-simulator/res/ccf")
        if not sim_separate_ccds:
            # Simulate xmm using one big ccd covering the fov.
            self.ccds.append(
                PnCcd(num=0, res_mult=res_mult, rotation=90.0, individual_ccd=False, ccf_path=ccf_path)
            )
        else:
            # Create the separate ccd.py's, the ccd.py alignment for higher resolutions is still off # TODO;
            print("WARNING: THE CCD ALIGNMENT FOR 2X AND 4X IS NOT CORRECT")
            quadrants_0_1 = [
                PnCcd(num=i, res_mult=res_mult, rotation=90.0, individual_ccd=True, ccf_path=ccf_path)
                for i in range(1, 7)
            ]
            quadrants_2_3 = [
                PnCcd(num=i, res_mult=res_mult, rotation=270.0, individual_ccd=True, ccf_path=ccf_path)
                for i in range(7, 13)
            ]
            self.ccds = quadrants_0_1.extend(quadrants_2_3)

    def save_xml(self, folder_location: Path):
        for ccd in self.ccds:
            ccd.save_xml(folder_location, self.base_name)

    def pprint_xml(self):
        for i, ccd in enumerate(self.ccds):
            print(f"----------- ccd{ccd.num:02d}.xml ---------------")
            ccd.pprint_xml()
            print("")
