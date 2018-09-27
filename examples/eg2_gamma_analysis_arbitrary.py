"""
gamma_analysis_arbitrary.py

Example of gamma analysis of distributions from abitrary sources

In this example, constructing instances of the Distribution class is
demonstrated. This allows any distribution to be imported that isn't explicitly
supported by the load_file() function.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import special
from flashgamma import Distribution, gamma_evaluation

# We need a 2D array of dose data. This array could come from any source we
# like. For this example, let's just create a simple square.
data = np.zeros((201, 201))
data[50:150, 50:150] = 2.0

# Now that we have the 2D dose data, we can create a reference Distribution
# object. If we set position=None, it will assume the origin is at the centre
# of the data.
ref_dist = Distribution(data, resolution=1, position=None)

# We also need an evaluated distribution. Let's just use the same square data.
eval_dist = Distribution(data, resolution=1, position=None)

# To make things interesting, let's give the reference distribution a
# small shift
ref_dist.translation = np.array([5, 5])

# We now need to reinterpolate the evaluated distribution at the new reference
# distribution positions. 
#
# For the DTA calculation, we also need to decide what resolution we want the 
# evaluated distribution to be. The greater the better (with diminishing 
# returns), but at the cost of increased calculation time. If it is too small, 
# accuracy will be lost due to the coarse discretisation. There is some debate 
# in the literature as to the optimal value; somewhere between 3 to 10 times
# finer than the size of the DTA criterion is probably acceptable. We will
# choose a DTA of 2 mm, so let's set the resolution to 3 points / mm (6 times
# finer than the DTA)
eval_dist = eval_dist.scale_grid(
    reference_distribution=ref_dist,
    new_resolution=3
)

# Now let's perform the gamma evaluation
gamma_dist, pass_rate = gamma_evaluation(
    ref_dist,
    eval_dist,
    delta_dose=2,
    delta_distance=2,
    threshold=0,
    local=False,
    pass_rate_only=False
)

# Finally, let's plot the results.
f1, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,4.2))
im1 = ax1.pcolormesh(ref_dist.position[0], ref_dist.position[1], ref_dist.data)
ax1.set_title('Reference Distribution')
ax1.set_xlabel('X (mm)')
ax1.set_ylabel('Y (mm)', labelpad=-5)
c1 = f1.colorbar(im1, ax=ax1)
c1.set_label('Dose (Gy)', rotation=270, labelpad=12)

im2 = ax2.pcolormesh(eval_dist.position[0], eval_dist.position[1], eval_dist.data)
ax2.set_title('Evaluated Distribution')
ax2.set_xlabel('X (mm)')
ax2.set_ylabel('Y (mm)', labelpad=-5)
c2 = f1.colorbar(im2, ax=ax2)
c2.set_label('Dose (Gy)', rotation=270, labelpad=12)

im3 = ax3.pcolormesh(gamma_dist.position[0], gamma_dist.position[1], gamma_dist.data, cmap='RdBu_r')
ax3.set_title('Gamma Evaluation (2%, 2mm): {:.1f}%'.format(pass_rate))
ax3.set_xlabel('X (mm)')
ax3.set_ylabel('Y (mm)', labelpad=-5)
c3 = f1.colorbar(im3, ax=ax3)
c3.set_label('Gamma Index', rotation=270)

plt.tight_layout()
plt.show()
