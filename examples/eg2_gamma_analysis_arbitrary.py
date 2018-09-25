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

# We also need an evaluated distribution. Let's just use the same data for now.
# However, for the DTA calculations to make sense, we need to increase its
# resolution. We can do this using .scale_grid()
eval_dist = Distribution(eval_data, resolution=2, position=None)

# To make thinhs interesting, let's give the the reference distribution a
# little twist
# ref_dist.rotation = np.pi / 8

ref_dist.translation = np.array([10, 10])
eval_dist = eval_dist.scale_grid(reference_distribution=ref_dist, scale=1)

# # Now let's have a look at the gamma evaluation
# gamma_dist = gamma_evaluation(
#     ref_dist,
#     eval_dist,
#     delta_dose=3,
#     delta_distance=3,
#     threshold=0,
#     local=False,
#     pass_rate_only=False
# )

# xx, yy = np.meshgrid(x, y)

# f = plt.figure()
# ax1 = f.gca(projection='3d')
# im1 = ax1.plot_surface(xx, yy, dist,cmap=cm.coolwarm)

# f = plt.figure()
# ax1 = f.gca()

# ax1.plot(eval_dist.position[1, :, 200], eval_dist.data[:, 200], label='eval')
# ax1.plot(ref_dist.position[1, :, 100], ref_dist.data[:, 100], label='ref')
# ax1.legend()

f = plt.figure()
ax1 = f.gca()
im1 = ax1.pcolormesh(ref_dist.position[0], ref_dist.position[1], ref_dist.data)#, cmap=cm.jet)
f.colorbar(im1, ax=ax1)

f = plt.figure()
ax1 = f.gca()
im1 = ax1.pcolormesh(eval_dist.position[0], eval_dist.position[1], eval_dist.data)#, cmap=cm.jet)
f.colorbar(im1, ax=ax1)

# f = plt.figure()
# ax1 = f.gca()
# im1 = ax1.pcolormesh(gamma_dist.position[0], gamma_dist.position[1], gamma_dist.data, cmap='RdBu_r')
# f.colorbar(im1, ax=ax1)

plt.show()
