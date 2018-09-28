import numpy as np
import csv
import math
from .distribution import Distribution


def load_file(filename, filetype):
    """Create a new Distribution from file data.

    Parses a file to create a new Distribution for analysis. Currently only
    ArcCheck measured files and SNC extracted files are fully supported.
    TODO(sam) Clean up film imports and do DICOM RT Dose importing.

    Parameters
    ----------
    filename : str
        Name of the file to be loaded.
    filetype : str
        Type of file to be loaded. Can be "arccheck" for ArcCheck measured
        files, "snc_extracted" for SNC extracted files, of "tiff" for .tiff
        images.

    Returns
    -------
    Distribution
        A new distrubution created with data values and locations parsed
        from the file information.
    """
    with open(filename, 'r') as f:

        if filetype == "arccheck":
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)[278:321]  # Grabs "Dose Counts" table
            data = np.array(rows)
            y = data[:-2:2, 0][::-1]  # Grab non-zero row/cols only
            x = data[-1, 2::2]
            xx, yy = np.meshgrid(x, y)
            position = np.stack([xx, yy])
            position = position.astype(float) * 10  # Convert cm to mm
            data = np.flipud(data[:-2:2, 2::2])  # Slice'n'dice to remove zeros
            data = data.astype(float)
            new_distribution = Distribution(
                data, position=position, resolution=0.1
            )

        elif filetype == "snc_extracted":
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)[5:]  # Snip off header
            data = np.array(rows)
            data[0, 0] = 0  # Clear Y\X table label
            data = data[:, :-1]
            data = data.astype(float)
            y = data[1:, 0][::-1]
            x = data[0, 1:]
            xx, yy = np.meshgrid(x, y)
            position = np.stack([xx, yy])
            data = np.flipud(data[1:, 1:])  # Slice'n'dice to remove labels
            new_distribution = Distribution(
                data, position=position, resolution=1
            )

        elif filetype == "dicom":
            raise NotImplementedError("DICOM import not implemented yet.")

        elif filetype == "tiff":
            # data = np.flipud(imread(filename)[:, :, 0]) # Take R channel only
            # # TODO: Fix hardcoded film calibration
            # data = np.log10(48000) - np.log10(data)
            # data = data.clip(min=0)
            # a = 12806.9125486
            # b = 2.9355787362
            # c = 3098.623464
            # d = -8.22823607708
            # data = a * np.power(data, b) + c * data + d
            # new_distribution = Distribution(data, resolution=2.95)
            raise NotImplementedError(".tiff import not implemented yet.")

        else:
            raise NotImplementedError("This file type is not yet supported.")

        return new_distribution

def create_distance_kernel(dta, resolution):
    """Create a distance kernel for gamma evaluations.

    Constructs a distance kernel for aiding calculations of the distance
    criteria component of gamma evaluations. This kernel only needs to be
    computed once for each distance to agreement criterion.

    Parameters
    ----------
    dta : float
        Distance to agreement crtierion, in mm.
    resolution : float
        Resolution of the evaluated distribution, in points / mm.

    Returns
    -------
    ndarray
        Returns a kernel to be fed into gamma evaluations.
    """
    eval_grids = math.ceil(dta * resolution)
    x = np.linspace(-eval_grids, eval_grids, 2 * eval_grids + 1)
    kernel = np.stack(np.meshgrid(x, x))
    kernel = np.abs(kernel)
    kernel /= resolution
    kernel = np.sum(kernel, axis=0)
    kernel = kernel**2
    kernel[np.sqrt(kernel) > dta] = 1000000000.0
    kernel = kernel / dta**2
    return kernel
