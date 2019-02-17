"""
gamma_surface.py

Example of a gamma surface analysis

In this example, the gamma pass rate is calculated for a large number of dose
difference and distance to agreement criteria. The results are presented in a
2D gamma surface plot.
"""
from os.path import dirname, abspath, join
import numpy as np
import matplotlib.pyplot as plt
from flashgamma import load_file, gamma, gamma_pass_rate

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
# For the DTA calculation, we also need to decide what resolution we want the
# evaluated distribution to be. If it is too poor, accuracy will be lost due to
# the coarse discretisation. There is some debate in the literature as to the
# optimal value; somewhere between 3 to 10 times finer than the size of
# the DTA criterion is probably acceptable. We will choose a smallest DTA of
# 1 mm, so let's set the resolution to 3 points / mm (3 times finer than the
# smallest DTA)
eval_dist = eval_dist.scale_grid(
    reference_distribution=ref_dist, new_resolution=3
)

# Define our dose difference and DTA criteria
dose_criteria = np.linspace(1.0, 3.0, 21)  # %
distance_criteria = np.linspace(1.0, 3.0, 21)  # mm

# Perform gamma evaluations in global mode, and no low dose threshold.
pass_rates = gamma_pass_rate(
    ref_dist,
    eval_dist,
    delta_dose=dose_criteria,
    delta_distance=distance_criteria,
    threshold=0,
    local=False
)

# Now let's plot the gamma surface.
f = plt.figure()
ax1 = plt.gca()

xx, yy = np.meshgrid(distance_criteria, dose_criteria, )
im1 = ax1.pcolormesh(xx, yy, pass_rates, shading='gouraud')

co1 = ax1.contour(distance_criteria, dose_criteria, pass_rates, colors='k', linewidth=1)
zc = co1.collections
plt.setp(zc, linewidth=1)
plt.clabel(co1, fontsize=10, fmt='%1.0f')

co2 = ax1.contour(distance_criteria, dose_criteria, pass_rates, colors='#f70c0c', levels=[95])
zc = co2.collections
plt.setp(zc, linewidth=2)
plt.clabel(co2, fontsize=12, fmt='%1.0f')

ax1.set_title('Gamma Pass Rates')
ax1.set_xlabel('Distance to agreement (mm)')
ax1.set_ylabel('Dose difference (%)')
c1 = plt.colorbar(im1, ax=ax1)
c1.set_label('Pass rate (%)', rotation=270, labelpad=12)

plt.show()