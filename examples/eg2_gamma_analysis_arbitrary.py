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
# like. For this example, let's create a 10 x 10 cm square that looks like it
# could have come from a linac beam. The function used here comes from !REF

# Define the necessary constants
eta = 1
T = 0.01
A = 0.173
B1 = 0.456
B2 = 2.892
x0 = 0
C = 0

# Create an array of 100 points and stick them into the function
x = np.linspace(-50, 50, 100)
D1 = (
    eta * ( T + (1 - T)*(A * special.erf(B1 * (x0 - x)) + (1 - A) *
    special.erf(B2 * (x0 - x)))) + C
)
D1 = D1 / 2 + 0.5

# Make a mirror image of the curve and stick the two ends together
D2 = D1[::-1]
D = np.append(D2, 1.0)
D = np.append(D, D1)

# Stack copies of the curve next to each along the x and y directions
data_x = D
for i in range(200):
    data_x = np.vstack((data_x, D))
data_y = D[np.newaxis].T
for i in range(200):
    data_y = np.hstack((data_y, D[np.newaxis].T))

# Combine the two directional stacks to get the final 2D surface
ref_data = data_x * data_y

# Do it all again at higher resolution for the evaluated data
x = np.linspace(-50, 50, 200)
D1 = (
    eta * ( T + (1 - T)*(A * special.erf(B1 * (x0 - x)) + (1 - A) *
    special.erf(B2 * (x0 - x)))) + C
)
D1 = D1 / 2 + 0.5
D2 = D1[::-1]
D = np.append(D2, 1.0)
D = np.append(D, D1)
data_x = D
for i in range(400):
    data_x = np.vstack((data_x, D))
data_y = D[np.newaxis].T
for i in range(400):
    data_y = np.hstack((data_y, D[np.newaxis].T))
eval_data = data_x * data_y


# Now that we have the 2D dose data, we can create a reference Distribution
# object. If we set position=None, it will assume the origin is at the centre
# of the data.
ref_dist = Distribution(ref_data, resolution=1, position=None)

# We also need an evaluated distribution. Let's use the higher resolution data
# we prepared earlier.
eval_dist = Distribution(eval_data, resolution=2, position=None)

# To make things interesting, let's give the reference distribution a
# small shift
ref_dist.translation = np.array([10, 10])

# We need to reinterpolate the evaluated distribution at the new reference
# distribution positions. 
#
# For the DTA calculation, we also need to decide what resolution we want the 
# evaluated distribution to be. The greater the better, but at the cost of 
# increased calculation time. If it is too small, accuracy will be lost due to 
# the coarse discretisation. There is some debate in the literature as to the 
# optimal value; somewhere between 3 to 10 times finer than the size of
# the DTA criterion is probably acceptable. We will choose a DTA of 2 mm, so
# let's set the resolution to 3 points / mm (6 times finer than the DTA)
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

# Finally let's plot the results.
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
