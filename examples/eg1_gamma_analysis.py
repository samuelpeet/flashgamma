"""
gamma_analysis.py

Example of a standard gamma analysis

In this example, a standard gamma analysis is performed on data acquired with
an arccheck. The basic workflow is demonstrated and the full gamma
distribution is plotted.
"""
from os.path import dirname, abspath, join
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from flashgamma import load_file, gamma

# Load the reference distribution. This is typically from a measurement, such
# as with an arccheck.
filename = "sample_arccheck_measurement.txt"
ref_dist = load_file(filename, "arccheck")

# Now load the evaluated distribution. This is typically plan data from a TPS.
# In this case we will load the plan data matching the arccheck diode positions
# extracted by SNC Patient.
filename = "sample_snc_extracted_data.snc"
eval_dist = load_file(filename, "snc_extracted")

# The two distributions have points at different positions and resolutions.
# We need to interpolate the evaluated distribution at the points of the
# reference distribution.
#
# For the DTA calculation, we need to decide what resolution we want the
# evaluated distribution to be. If it is too poor, accuracy will be lost due to
# the coarse discretisation. There is some debate in the literature as to the
# optimal value; somewhere between 3 to 10 times finer than the size of
# the DTA criterion is probably acceptable. We will choose a DTA of 2 mm, so
# let's set the resolution to 3 points / mm (6 times finer than the DTA)
eval_dist = eval_dist.scale_grid(
    reference_distribution=ref_dist, new_resolution=3
)

# Perform a gamma evaluation in global mode, with 2% dose difference, 2 mm DTA,
# and no low dose threshold.
gamma_dist, pass_rate = gamma(
    ref_dist,
    eval_dist,
    delta_dose=2,
    delta_distance=2,
    threshold=0,
    local=False
)

# Finally let's plot the results.
f = plt.figure(figsize=(14,7))
gs = gridspec.GridSpec(2, 4)
gs.update(
    left=0.05,
    right=1.00,
    top=0.95,
    bottom=0.08,
    wspace=0.11,
    hspace=0.40
)
ax1 = plt.subplot(gs[0, :2], )
ax2 = plt.subplot(gs[0, 2:])
ax3 = plt.subplot(gs[1, 1:3])

im1 = ax1.pcolormesh(ref_dist.position[0], ref_dist.position[1], ref_dist.data)
ax1.set_title('Reference Distribution')
ax1.set_xlabel('X (mm)')
ax1.set_ylabel('Y (mm)', labelpad=-5)
c1 = plt.colorbar(im1, ax=ax1)
c1.set_label('Dose (cGy)', rotation=270, labelpad=12)

im2 = ax2.pcolormesh(eval_dist.position[0], eval_dist.position[1], eval_dist.data)
ax2.set_title('Evaluated Distribution')
ax2.set_xlabel('X (mm)')
ax2.set_ylabel('Y (mm)', labelpad=-5)
c2 = plt.colorbar(im2, ax=ax2)
c2.set_label('Dose (cGy)', rotation=270, labelpad=12)

gamma_dist._data[gamma_dist.data == np.inf] = 0
cold = gamma_dist.data < -1
cold_x = gamma_dist.position[0, cold]
cold_y = gamma_dist.position[1, cold]
hot = gamma_dist.data > 1
hot_x = gamma_dist.position[0, hot]
hot_y = gamma_dist.position[1, hot]
im3 = ax3.pcolormesh(eval_dist.position[0], eval_dist.position[1], eval_dist.data, alpha=0.1, cmap='Greys')
ax3.plot(cold_x, cold_y, 'vb', markersize=5)
ax3.plot(hot_x, hot_y, '^r', markersize=5)
ax3.set_title('Gamma Evaluation (2%, 2mm): {:.1f}%'.format(pass_rate))
ax3.set_xlabel('X (mm)')
ax3.set_ylabel('Y (mm)', labelpad=-5)

plt.show()