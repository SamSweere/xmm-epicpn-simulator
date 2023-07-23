RESOLUTION = 402
PIXEL_SIZE = 150e-4  # based of: https://www.cosmos.esa.int/web/xmm-newton/boundaries-pn
# cdelt give the pixel sizes in degrees
# cdelt from XMM_MISCDATA_0022.CCF PLATE_SCALE_X, the unit is in arsec, arsec to degree by deciding it by 3600
# TODO: this does not cover the whole xmm background
CDELT = (4.12838 / 3600)
