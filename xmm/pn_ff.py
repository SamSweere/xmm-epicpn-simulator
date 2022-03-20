from utils.log import plog
from xmm.pn_ccd import XmmPnCcd


class XmmPnFf:
    def __init__(self, res_mult=1, noise=True, sim_separate_ccds=False, verbose=True):
        self.verbose = verbose
        self.base_name = "fullframe_thinfilter"

        # Based on EPN_LINCOORD_0009.CCF
        # The location of the center of each ccd
        # self.xyrval = {
        #     '01': (-4.832e-03, 15.1455e-03),
        #     '02': (-14.626e-03, 15.1455e-03),
        #     '03': (-24.42e-03, 15.1455e-03),
        #     '04': (4.982e-03, 15.1455e-03),
        #     '05': (14.776e-03, 15.1455e-03),
        #     '06': (24.57e-03, 15.1455e-03),
        #     '07': (4.832e-03, -15.1455e-03),
        #     '08': (14.626e-03, -15.1455e-03),
        #     '09': (24.42e-03, -15.1455e-03),
        #     '10': (-4.982e-03, -15.1455e-03),
        #     '11': (-14.776e-03, -15.1455e-03),
        #     '12': (-24.57e-03, -15.1455e-03),
        # }

        self.xyrval = {
            '01': (15.1455e-03, -4.832e-03),
            '02': (15.1455e-03, -14.626e-03),
            '03': (15.1455e-03, -24.42e-03),
            '04': (15.1455e-03, 4.982e-03),
            '05': (15.1455e-03, 14.776e-03),
            '06': (15.1455e-03, 24.57e-03),
            '07': (-15.1455e-03, 4.832e-03),
            '08': (-15.1455e-03, 14.626e-03),
            '09': (-15.1455e-03, 24.42e-03),
            '10': (-15.1455e-03, -4.982e-03),
            '11': (-15.1455e-03, -14.776e-03),
            '12': (-15.1455e-03, -24.57e-03),
        }

        # Optical offset based on XMM_MISCDATA_0022.CCF:
        OPTICS_X = 23
        OPTICS_Y = 183
        # OPTICS_CCD = 4
        CCD_DEFAULT_PIXEL_WIDTH = 64
        CCD_DEFAULT_PIXEL_HEIGHT = 200
        CCD_DEFAULT_PIXEL_SIZE = 0.15e-03

        DEFAULT_IMG_WIDTH = 403
        DEFAULT_IMG_HEIGHT = 411

        self.x_optical_offset = -1.0 * (self.xyrval['04'][0] + (
                    CCD_DEFAULT_PIXEL_HEIGHT / 2.0 - OPTICS_Y) * CCD_DEFAULT_PIXEL_SIZE)  # 0.0026955
        self.y_optical_offset = -1.0 * (self.xyrval['04'][1] - (
                    OPTICS_X - CCD_DEFAULT_PIXEL_WIDTH / 2.0) * CCD_DEFAULT_PIXEL_SIZE)  # 0.003632

        # Calculate the ccd pixel width and height and also the pixel size
        if sim_separate_ccds:
            ccd_pixel_width = CCD_DEFAULT_PIXEL_WIDTH * res_mult
            ccd_pixel_height = CCD_DEFAULT_PIXEL_HEIGHT * res_mult
        else:
            ccd_pixel_width = DEFAULT_IMG_WIDTH * res_mult
            ccd_pixel_height = DEFAULT_IMG_HEIGHT * res_mult
        p_delt = CCD_DEFAULT_PIXEL_SIZE / float(res_mult)  # The size of one pixel

        # cdelt give the pixel sizes in degrees
        # cdelt from XMM_MISCDATA_0022.CCF PLATE_SCALE_X, the unit is in arsec, arsec to degree by deciding it by 3600
        # The cdelt is scaled in order to have the same layout as the real images. I did not find the reason why this is needed
        # The scaling factor is determined by the pixel difference on the x-axis
        self.cdelt_x = 4.0 / 3600 / float(res_mult)
        self.cdelt_y = 4.0 / 3600 / float(res_mult)

        # focal_length from from XMM_MISCDATA_0022.CCF XRT3, from mm to meters
        focal_length = 7493.2 / 1000  # meter

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

        # Old method, resolution invariant but is only correct for 1x scaling
        # crpix1_corr = -4*res_mult
        # crpix2_corr = 3*res_mult

        # self.crpix1 = self.total_pixel_height/2.0 - self.y_optical_offset / p_delt + crpix1_corr
        # self.crpix2 = self.total_pixel_width/2.0 - self.x_optical_offset/p_delt + crpix2_corr

        # The fov is bigger than the determined by the biggest axis and the fov per pixel
        # based of XMM_MISCDATA_0022.CCF XRT3 FOV_RADIUS (0.5). The reason it is bigger is to have some margin
        # the 0.5 fov created too small images
        fov = 0.8

        # Readout time
        readout_wait_time = 68.75e-3

        # Setting this to 0.0 eliminates out of time events
        lineshift_wait_time = 0.0  # 23.04e-6 # /res_mult

        # Create the global settings dictionary to be passed to the individual ccd's
        self.g_settings = {
            'res_mult': res_mult,
            'base_name': self.base_name,
            'p_delt': p_delt,
            'p_delt_default': CCD_DEFAULT_PIXEL_SIZE,
            'fov': fov,
            'focal_length': focal_length,
            'width': ccd_pixel_width,
            'height': ccd_pixel_height,
            'readout_wait_time': readout_wait_time,
            'lineshift_wait_time': lineshift_wait_time,
            'noise': noise
        }

        self.ccds = []

        if not sim_separate_ccds:
            # Simulate xmm using one big ccd covering the fov.
            ccd = XmmPnCcd(num=0, g_settings=self.g_settings, xrval=self.x_optical_offset,
                           yrval=self.y_optical_offset, rotation=90.0, verbose=verbose)

            self.ccds.append(ccd)
        else:
            # Create the separate ccd's, the ccd alignment for higher resolutions is still off # TODO;
            print("WARNING: THE CCD ALIGNMENT FOR 2X AND 4X IS NOT CORRECT")
            self.create_ccds()

    def create_ccds(self):
        # Note the default settings are for ccd01
        # Quadrant 0
        self.ccd01 = XmmPnCcd(num=1, g_settings=self.g_settings, xrval=self.xyrval['01'][0] + self.x_optical_offset,
                              yrval=self.xyrval['01'][1] + self.y_optical_offset, rotation=90.0, verbose=self.verbose)
        self.ccd02 = XmmPnCcd(num=2, g_settings=self.g_settings, xrval=self.xyrval['02'][0] + self.x_optical_offset,
                              yrval=self.xyrval['02'][1] + self.y_optical_offset, rotation=90.0, verbose=self.verbose)
        self.ccd03 = XmmPnCcd(num=3, g_settings=self.g_settings, xrval=self.xyrval['03'][0] + self.x_optical_offset,
                              yrval=self.xyrval['03'][1] + self.y_optical_offset, rotation=90.0, verbose=self.verbose)

        # Quadrant 1
        self.ccd04 = XmmPnCcd(num=4, g_settings=self.g_settings, xrval=self.xyrval['04'][0] + self.x_optical_offset,
                              # - self.g_settings['p_delt'], # (self.g_settings['res_mult'] - 1) * self.g_settings['p_delt'], #TODO
                              yrval=self.xyrval['04'][1] + self.y_optical_offset, rotation=90.0, verbose=self.verbose)
        self.ccd05 = XmmPnCcd(num=5, g_settings=self.g_settings, xrval=self.xyrval['05'][0] + self.x_optical_offset,
                              # - self.g_settings['p_delt'],# (self.g_settings['res_mult'] - 1) * self.g_settings['p_delt'], #TODO
                              yrval=self.xyrval['05'][1] + self.y_optical_offset, rotation=90.0, verbose=self.verbose)
        self.ccd06 = XmmPnCcd(num=6, g_settings=self.g_settings, xrval=self.xyrval['06'][0] + self.x_optical_offset,
                              # - self.g_settings['p_delt'],# (self.g_settings['res_mult'] - 1) * self.g_settings['p_delt'], #TODO
                              yrval=self.xyrval['06'][1] + self.y_optical_offset, rotation=90.0, verbose=self.verbose)

        # Quadrant 2
        self.ccd07 = XmmPnCcd(num=7, g_settings=self.g_settings,
                              xrval=self.xyrval['07'][0] + self.x_optical_offset + 2 * self.g_settings[
                                  'p_delt_default'],
                              yrval=self.xyrval['07'][1] + self.y_optical_offset, rotation=270.0,
                              verbose=self.verbose)  # + 2*self.g_settings['p_delt_default']
        self.ccd08 = XmmPnCcd(num=8, g_settings=self.g_settings,
                              xrval=self.xyrval['08'][0] + self.x_optical_offset + 2 * self.g_settings[
                                  'p_delt_default'],
                              yrval=self.xyrval['08'][1] + self.y_optical_offset, rotation=270.0,
                              verbose=self.verbose)  # + 2*self.g_settings['p_delt']
        self.ccd09 = XmmPnCcd(num=9, g_settings=self.g_settings,
                              xrval=self.xyrval['09'][0] + self.x_optical_offset + 2 * self.g_settings[
                                  'p_delt_default'],
                              yrval=self.xyrval['09'][1] + self.y_optical_offset, rotation=270.0,
                              verbose=self.verbose)  # + 2*self.g_settings['p_delt']

        # Quadrant 3
        self.ccd10 = XmmPnCcd(num=10, g_settings=self.g_settings,
                              xrval=self.xyrval['10'][0] + self.x_optical_offset + 2 * self.g_settings[
                                  'p_delt_default'],
                              yrval=self.xyrval['10'][1] + self.y_optical_offset, rotation=270.0,
                              verbose=self.verbose)  # 0.5*self.g_settings['p_delt_default'] seems to work for 2x res_mult
        self.ccd11 = XmmPnCcd(num=11, g_settings=self.g_settings,
                              xrval=self.xyrval['11'][0] + self.x_optical_offset + 2 * self.g_settings[
                                  'p_delt_default'],
                              yrval=self.xyrval['11'][1] + self.y_optical_offset, rotation=270.0,
                              verbose=self.verbose)  # + 2*self.g_settings['p_delt']
        self.ccd12 = XmmPnCcd(num=12, g_settings=self.g_settings,
                              xrval=self.xyrval['12'][0] + self.x_optical_offset + 2 * self.g_settings[
                                  'p_delt_default'],
                              yrval=self.xyrval['12'][1] + self.y_optical_offset, rotation=270.0,
                              verbose=self.verbose)  # + 2*self.g_settings['p_delt']

        self.ccds = [self.ccd01, self.ccd02, self.ccd03, self.ccd04, self.ccd05, self.ccd06, self.ccd07, self.ccd08,
                     self.ccd09, self.ccd10, self.ccd11, self.ccd12]

    def save_xml(self, folder_location):
        for i in range(len(self.ccds)):
            ccd = self.ccds[i]
            ccd.save_xml(folder_location)
        plog("All ccd's saved", verbose=self.verbose)

    def pprint_xml(self):
        for i in range(len(self.ccds)):
            ccd = self.ccds[i]
            print("----------- ccd{:02d}".format(ccd.num) + ".xml ---------------")
            ccd.pprint_xml()
            print("")
