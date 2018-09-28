import numpy as np
import math
from .distribution import Distribution
from .helpers import create_distance_kernel
from .gamma_evaluation_c import gamma_evaluation_c

def difference_between(dist1, dist2, kind="relative"):
    """Compute the difference between two distributions.

    Calculate either the absolute or relative difference between the values
    in the two distributions. Both distributions must have matching grid sizes
    and positions.

    Parameters
    ----------
    distribution1 : Distribution
        The first distribution to calculate the difference between.
    distribution2 : Distribution
        The second distribution to calculate the difference between.
    kind : str
        The type of comparison to be made. Either "relative" or "absolute".

    Returns
    -------
    Distribution
        A new distribution with data values corresponding to the difference
        between the two input distributions.
    """
    if kind == "absolute":
        new_data = dist2.data - dist1.data
    elif kind == "relative":
        new_data = np.copy(dist2.data)
        new_data = (
            (dist1.data - dist2.data) /
            np.max(dist2.data) * 100
        )
        new_data = np.nan_to_num(new_data)

    new_distribution = Distribution(
        new_data,
        resolution=dist2.resolution,
        position=dist2.position
    )
    return new_distribution

def gamma(r_dist, e_dist, delta_dose=3, delta_distance=3,
                        threshold=0, local=False):
    """Perform a gamma evaluation between a reference distribution and an
    evaluated distribution.

    Performs a gamma analysis between a reference distribution and an
    evaluated distribution. Parts of this code are inspired by Christopher
    Poole's pygamma sofware at https://github.com/christopherpoole/pygamma

    Parameters
    ----------
    r_dist : Distribution
        The reference distribution. Typically this would be measured dose
        data.
    e_dist : Distribution
        The distribution to evaluate gamma indices against. Typically this
        would be radiotherapy plan data.
    delta_dose : float
        Dose difference criterion, in %.
    delta_distance : float
        Distance to agreement criterion, in mm.
    threshold : float
        Fraction of the maximum evaluated distribution value under which
        gamma calculations will be skipped, in %.
    local : bool
        Perform a local gamma evaluation. If false, perform a Van Dyk global
        evaluation instead.

    Returns
    -------
    Distribution
        A new distrubution with data values corresponding to the calculated
        gamma index at each position.
    float
        The overall pass rate of the gamma evaluation, in %.
    """
    # Create gamma map and set all points to a passing index
    gamma_map = np.ones_like(r_dist.data)
    gamma_map[:, :] = np.inf

    # Determine how many grid points wide the dta search needs to be
    eval_grids = math.ceil(delta_distance * e_dist.resolution)

    # Construct the dta kernel
    x = np.linspace(-eval_grids, eval_grids, 2 * eval_grids + 1)
    kernel = np.stack(np.meshgrid(x, x))
    kernel = np.abs(kernel)
    kernel /= e_dist.resolution
    kernel = np.sum(kernel, axis=0)
    kernel = kernel**2
    kernel[np.sqrt(kernel) > delta_distance] = np.inf
    kernel = kernel / delta_distance**2

    # Iterate through reference distribution and compute gamma index for
    # each point
    # TODO(sam) Deal with boundary conditions better
    ref_grids = math.ceil(delta_distance * r_dist.resolution)
    max_e_dose = np.max(e_dist.data)
    for i in range(ref_grids, r_dist.data.shape[0] - ref_grids):
        for j in range(ref_grids, r_dist.data.shape[1] - ref_grids):

            # Dose at reference grid point
            D_r = r_dist.data[i, j]

            # Check if under threshold
            if D_r < max_e_dose * threshold / 100:
                gamma_map[i, j] = np.inf
                continue

            # Find correspoding point in evaluated distribution
            res_scale = int(e_dist.resolution / r_dist.resolution)
            c_r, c_c = (i * res_scale, j * res_scale)

            # Slice evaluated distribution to same shape as distance kernel
            min_r, max_r = (c_r - eval_grids, c_r + eval_grids)
            min_c, max_c = (c_c - eval_grids, c_c + eval_grids)
            e_slice = e_dist.data[min_r:max_r+1, min_c:max_c+1]
            assert e_slice.shape == kernel.shape, \
                "e_slice and kernel must be the same shape"

            # Difference between doses at evaluated points and reference point
            e_slice = e_slice - D_r

            # Handle local or global gamma criteria
            if local:
                denom = D_r * delta_dose / 100
            else:
                denom = max_e_dose * delta_dose / 100

            # Compute dose half of gamma value
            e_slice = e_slice**2 / denom**2

            # Apply kernel (distace half of gamma value)
            e_slice = e_slice + kernel

            # Minimise result to get square of gamma index
            gamma_index = np.min(e_slice)

            # Store gamma index
            gamma_map[i, j] = np.sqrt(gamma_index)

            # Encode information as to whether point was hot or cold
            if D_r < e_dist.data[c_r, c_c]:
                gamma_map[i, j] = gamma_map[i, j] * -1  # Point is cold

    num_passing = np.sum(np.abs(gamma_map) <= 1.0)
    total_points = gamma_map.size - np.sum(gamma_map == np.inf)
    pass_rate = num_passing / total_points * 100

    new_distribution = Distribution(
        gamma_map,
        resolution=r_dist.resolution,
        position=r_dist.position
    )
    return new_distribution, pass_rate

def gamma_pass_rate(r_dist, e_dist, delta_dose=3, delta_distance=3,
                        threshold=0, local=False):
    """Perform a gamma evaluation between a reference distribution and an
    evaluated distribution, and return the pass rate only

    Performs a gamma analysis between a reference distribution and an
    evaluated distribution. Parts of this code are inspired by Christopher
    Poole's pygamma sofware at https://github.com/christopherpoole/pygamma.
    Returns a pass rate only, allowing a more efficient gamma analysis
    algorithm implementation.

    Parameters
    ----------
    r_dist : Distribution
        The reference distribution. Typically this would be measured dose
        data.
    e_dist : Distribution
        The distribution to evaluate gamma indices against. Typically this
        would be radiotherapy plan data.
    delta_dose : float
        Dose difference criterion, in %. May be a single value or an ndarray.
    delta_distance : float
        Distance to agreement criterion, in mm. May be a single value or an
        ndarray.
    threshold : float
        Fraction of the maximum evaluated distribution value under which
        gamma calculations will be skipped, in %.
    local : bool
        Perform a local gamma evaluation. If false, perform a Van Dyk global
        evaluation instead.

    Returns
    -------
    float
        ndarray of pass rates, in %.
    """
    if isinstance(delta_dose, np.ndarray):
        pass
    elif isinstance(float(delta_dose), float):
        delta_dose = np.array([delta_dose])
    else:
        assert False, \
        "Dose criteria must be a single float or ndarray of floats"

    if isinstance(delta_distance, np.ndarray):
        pass
    elif isinstance(float(delta_distance), float):
        delta_distance = np.array([delta_distance])
    else:
        assert False, \
        "DTA criteria must be a single float or ndarray of floats"

    local = 1 if local else 0
    e_dist_max = np.max(e_dist.data)
    pass_rates = np.zeros((len(delta_dose), len(delta_distance)))

    for xi, xv in enumerate(delta_distance):

        # Construct kernel only once per DTA criterion
        kernel = create_distance_kernel(delta_distance[xi], e_dist.resolution)

        for yi, yv in enumerate(delta_dose):

            # Perform gamma anlysis
            pass_rates[yi, xi] = gamma_evaluation_c(
                r_dist.data, r_dist.resolution,
                e_dist.data, e_dist.resolution, e_dist_max,
                np.ones_like(r_dist.data),
                kernel,
                np.ones_like(kernel),
                d_dose=float(delta_dose[yi]), d_dist=float(delta_distance[xi]),
                thresh=float(threshold), local=int(local)
            )
    return pass_rates

def maximum_allowed_dose_difference(r_dist, e_dist, delta_dose=3,
                                    delta_distance=3, threshold=0, simple=True,
                                    acceptance_region='box',
                                    pass_rate_only=False, normalised=False):
    """Perform a maximum allowed dose difference evaluation between a
    reference distribution and an evaluated distribution.

    Perform a maximum allowed dose difference analysis between a reference
    distribution and an evaluated distribution. Further explanation
    of this analysis and algorithm can be found in "Jiang et al. On dose
    distribution comparison. Phys. Med. Biol. 51 (2006) 759–776".

    Parameters
    ----------
    r_dist : Distribution
        The reference distribution. Typically this would be measured dose
        data.
    e_dist : Distribution
        The distribution to evaluate gamma indices against. Typically this
        would be radiotherapy plan data.
    delta_dose : float
        Dose difference criterion, in %.
    delta_distance : float
        Distance to agreement criterion, in mm.
    threshold : float
        Fraction of the maximum evaluated distribution value under which
        gamma calculations will be skipped, in %.
    simple : bool
        If true, linearly interpolate the evaluated distribution at each
        reference point based on the gradient (Eq. 37 - 39). If false, perform
        the full search (Eq. 33 - 35).
    acceptance_region : str
        Define the acceptance region shape for the analysis. Options are
        'composite', 'box', and 'gamma'.
    local : bool
        Perform a local gamma evaluation. If false, perform a Van Dyk global
        evaluation instead.
    pass_rate_only : bool
        Return the overall pass rate only. If false, return the entire MADD
        map as a new distribution.
    pass_fail_map : bool
        Return 1 or 0 for each data point corresponding to passing or failing
        the analysis at that point.
    normalised: bool
        If true, returns a normalised equivalent dose distribution (Eq. 22).

    Returns
    -------
    float
        If pass_rate_only is true, returns the overall pass rate of the
        gamma evalution, in %.
    Distribution
        A new distrubution with data values corresponding to the calculated
        gamma index at each position.
    """
    # delta_dose = delta_dose / 100
    x = e_dist.position[0, 0, :]
    y = e_dist.position[1, :, 0]

    if simple:
        grad_e_dist = np.gradient(e_dist.data, y, x)
        norm_grad_e_dist = np.sqrt(
            np.power(grad_e_dist[0], 2) +
            np.power(grad_e_dist[1], 2)
        )
        beta = norm_grad_e_dist * delta_distance / delta_dose
    else:
        # TODO(sam) Implement full MADD calculation
        raise NotImplementedError(
            "Full MADD calculation not implemented yet!"
        )

    if acceptance_region == 'composite':
        madd = delta_dose * np.maximum(np.ones(beta.shape), beta)
    elif acceptance_region == 'gamma':
        beta2 = np.power(beta, 2)
        madd = delta_dose * (1 + 2 * beta2) / np.sqrt(1 + 4 * beta2)
    elif acceptance_region == 'box':
        madd = delta_dose * (1 + beta)
    else:
        raise AssertionError(
            "Acceptance region must be 'composite', 'box', or 'gamma'"
        )

    result = madd

    if normalised:
        dd = r_dist.difference_between(e_dist, kind="absolute").data
        result = delta_dose * dd / madd  #NDD

    if pass_rate_only:
        # TODO: Implement pass_rate_only
        raise NotImplementedError(
            "pass_rate_only for MADD calculation not implemented yet!"
        )

    if threshold:
        # TODO: Implement MADD thresholding
        assert NotImplementedError(
            "MADD thresholding not implemented yet!"
        )

    new_distribution = Distribution(
        result,
        resolution=r_dist.resolution,
        position=r_dist.position
    )

    return new_distribution
