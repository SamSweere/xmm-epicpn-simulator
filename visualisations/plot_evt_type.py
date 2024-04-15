import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

# Load the combined eventsfile
hdu = fits.open("/home/sam/Documents/ESA/data/sim/tmp/grid_point_sources_1_mult_10_0ks/ccd_combined_evt.fits")
# types = np.sort(hdu['EVENTS'].data['TYPE'])
unique, counts = np.unique(hdu["EVENTS"].data["TYPE"], return_counts=True)
perc = counts / np.sum(counts)

plt.bar(unique, perc)

# plt.hist(, bins=np.arange(14)-0.5)
plt.xlabel("Pattern Type")
plt.ylabel("Percentage")
plt.xticks(unique)
plt.yticks(np.arange(0, 0.5, 0.05))
# plt.savefig("/home/sam/Documents/ESA/thesis/figures/combined_evt_pattern.pdf")

plt.show()
