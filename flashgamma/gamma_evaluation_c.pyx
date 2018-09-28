cimport cython
cimport numpy as npc
from libc.math cimport ceil, sqrt
from cython.view cimport array as cvarray

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def gamma_evaluation_c(npc.ndarray[npc.float64_t, ndim=2] r_data, double r_res,
                       npc.ndarray[npc.float64_t, ndim=2] e_data, double e_res,
                       double e_max_value,
                       npc.ndarray[npc.float64_t, ndim=2] gamma_map,
                       npc.ndarray[npc.float64_t, ndim=2] kernel,
                       npc.ndarray[npc.float64_t, ndim=2] e_slice,
                       double d_dose=3, double d_dist=3, double thresh=0,
                       int local=0):
    """Perform a gamma evaluation between two distributions.

    Performs a gamma analysis between a reference distribution and an
    evaluated distribution. To maximise the execution speed, all ndarray
    creation must be performed outside of this function and passed in as
    parameters.

    Parameters
    ----------
    r_data : ndarray
        The reference distribution data. Typically this would be measured data.
    r_res : double
        The resolution of the reference distribution, in points / mm.
    e_data : ndarray
        The evaluated distribution data. Typically this would be plan data.
    e_res : double
        The resolution of the evaluated distribution, in points / mm.
    e_max_value : double
        The maximum value in the evaluated distribution.
    gamma_map: ndarray
        ndarray of the same shape as the reference distribution data. It
        becomes filled with computed gamma indices as this function
        executes.
    kernel: ndarray
        TODO: Update this description
    e_slice: ndarray
        ndarray of the same shape as the distance kernel. This array becomes
        filled with a slice of the evaluated distribution data surrounding each
        reference point during execution.
    d_dose : float
        Dose difference criterion, in %.
    d_dist : float
        Distance to agreement criterion, in mm.
    thresh : float
        Fraction of the maximum evaluated distribution value under which
        gamma calculations will be skipped, in %.
    local : int
        Perform a local gamma evaluation. If false, perform a Van Dyk global
        evaluation instead.

    Returns
    -------
    float
        The overall pass rate of the gamma evalution, in %.
    """
    cdef unsigned int i, j, m, n, loop_count
    cdef unsigned int c_r, c_c
    cdef unsigned int min_r, max_r, min_c, max_c
    cdef unsigned int emax
    cdef int res_scale, r_grids, e_grids
    cdef double D_r, denom
    cdef double min_gamma = 10000000.0
    cdef int total_points
    cdef int num_passing = 0

    r_grids = int(ceil(d_dist * r_res))
    e_grids = int(ceil(d_dist * e_res))
    res_scale = int(e_res / r_res)
    emax = 2 * e_grids + 1

    # Iterate through ref distribution and compute gamma index at each point
    # TODO: Deal with boundary conditions better
    for i in range(r_grids, <int>r_data.shape[0] - r_grids):
        for j in range(r_grids, <int>r_data.shape[1] - r_grids):

            # Dose at reference grid point
            D_r = <double>r_data[i, j]

            # Check if under threshhold threshold
            if D_r < e_max_value * thresh / 100.0:
                gamma_map[i, j] = 999999.0
                continue

            # Handle local or global gamma criteria
            if local:
                denom = D_r * d_dose / 100
            else:
                denom = e_max_value * d_dose / 100

            # Find correspoding point in evaluated distribution
            c_r = i * res_scale
            c_c = j * res_scale

            # Perform calculations over relevent part of distribution
            min_r = c_r - e_grids
            min_c = c_c - e_grids
            for m in range(emax):
                for n in range(emax):

                    # Difference between doses at evaluated points and reference point
                    e_slice[m, n] = e_data[min_r + m, min_c + n] - D_r

                    # Compute dose half of gamma value
                    e_slice[m, n] = e_slice[m, n]**2 / denom**2

                    # Apply kernel (distace half of gamma value)
                    e_slice[m, n] = e_slice[m, n] + kernel[m, n]

                    #print(e_slice[m, n])
                    if e_slice[m, n] < min_gamma:
                        min_gamma = e_slice[m, n]

            # Store gamma index
            gamma_map[i, j] = sqrt(min_gamma)

            min_gamma = 10000000.0

    total_points = gamma_map.shape[0] * gamma_map.shape[1]
    for i in range(<int>gamma_map.shape[0]):
        for j in range(<int>gamma_map.shape[1]):
            if gamma_map[i, j] <= 1.0:
                num_passing = num_passing + 1
            elif gamma_map[i, j] == 999999.0:
                total_points = total_points - 1

    return <double>num_passing / <double>total_points * 100.0