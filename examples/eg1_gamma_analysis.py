"""
gamma_analysis.py

Example of a standard gamma analysis

In this example, a standard gamma analysis is performed on data acquired with
an arccheck. The basic workflow is demonstrated and the full gamma
distribution is plotted.
"""
import matplotlib.pyplot as plt
from flashgamma import load_file, gamma_evaluation

# Load the reference distribution. This is typically from a measurement, such
# as with an arccheck
filename = "sample_arccheck_measurement.txt"
ref_dist = load_file(filename, "arccheck")

# Load the evaluated distribution. This is typically plan data from a TPS. In
# this case we will load the plan data matching the arccheck diode positions
# extracted by SNC Patient.
filename = "sample_snc_extracted_data.snc"
eval_dist = load_file(filename, "snc_extracted")

# These two distributions have points at different positions and resolutions.
# We need to interpolate the evaluated distribution at the points of the
# reference distribution. 
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

# Perform a gamma evaluation in global mode, with 2% dose difference, 2 mm DTA,
# and no low dose threshold. This will return a full gamma distribution. If we
# only wanted the pass rate, we would set pass_rate_only=True.
gamma_dist, pass_rate = gamma_evaluation(
    ref_dist,
    eval_dist,
    delta_dose=2,
    delta_distance=2,
    threshold=0,
    local=False,
    pass_rate_only=False
)