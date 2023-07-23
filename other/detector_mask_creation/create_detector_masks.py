import numpy as np
from astropy.io import fits

with fits.open('pn_expo_500_2000_detxy.ds') as hdu:
    # now make it a mask image
    #
    mask = hdu[0].data > 0.0
    #
    hdu_out = hdu.copy()
    hdu_out[0].data = mask.astype(np.ubyte)
    hdu_out[0].scale('ubyte')

    hdu_out.writeto('pn_mask_500_2000_detxy_1x.ds', overwrite=True)

    # Upscale the mask image 2x
    hdu_out[0].data = hdu_out[0].data.repeat(2, axis=0).repeat(2, axis=1)
    hdu_out.writeto('pn_mask_500_2000_detxy_2x.ds', overwrite=True)

    # Upscale the mask image 4x
    hdu_out[0].data = hdu_out[0].data.repeat(2, axis=0).repeat(2, axis=1)
    hdu_out.writeto('pn_mask_500_2000_detxy_4x.ds', overwrite=True)
