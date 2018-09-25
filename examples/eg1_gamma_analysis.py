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
# reference distribution. For the DTA calculation, we also need to decide how
# much greater we want the resolution of the evaluated distribution than the
# reference distribution. The greater the better, but at the cost of increased
# calculation time. If it is too small, accuracy will be lost due to the coarse
# discretisation. There is some debate in the literature as to the optimal
# value; somewhere between 3 and 10 is probably acceptable. Here we choose 5.
eval_dist = eval_dist.scale_grid(
    reference_distribution=ref_dist,
    scale=5
)

# Perform a gamma evaluation in global mode, with 2% dose difference, 2 mm DTA,
# and no low dose threshold. This will return a full gamma distribution. If we
# only wanted the pass rate, we would set pass_rate_only=True.
gamma_dist = gamma_evaluation(
    ref_dist,
    eval_dist,
    delta_dose=2,
    delta_distance=2,
    threshold=0,
    local=False,
    pass_rate_only=False
)