import copy
import numpy as np
from scipy import interpolate


class Distribution(object):
    """ Generic representation of an aribtrary distribution.

    This class implements a generic representation of a distribution. Its
    purpose is to faciliate comparisons between distributions that may come
    from different sources and have different resolutions.

    Parameters
    ----------
    data : ndarray
        Array of floats corresponding to the values of the distribution at each
        point.
    resolution : float
        Resolution of the dataset, in units of points / mm.
    position : ndarray
        Spatial position of each point in the dataset.
    """

    def __init__(self, data, resolution=1, position=None):
        self._data = data  # Point values as numpy array
        self._resolution = resolution  # Points / mm
        self._translation = np.array([0, 0])  # dx, dy
        self._rotation = 0  # Radians

        if position is None:
            # Create array containing position of each point in space. It is
            # assumed that the centre of the distribution is located at the
            # origin.
            x = range(self._data.shape[1])
            y = range(self._data.shape[0])
            xx, yy = np.meshgrid(x, y)
            self._position = np.stack([xx, yy])
            self._position = self._position.astype(float)
            # Convert from points to mm
            self._position = self._position / resolution
            # Move origin to centre of dataset
            x_offset = (self._position[0, 0, -1] - self._position[0, 0, 0]) / 2
            self._position[0, :, :] = self._position[0, :, :] - x_offset
            y_offset = (self._position[1, -1, 0] - self._position[1, 0, 0]) / 2
            self._position[1, :, :] = self._position[1, :, :] - y_offset
        else:
            self._position = position

    @property
    def data(self):
        """ndarray: Dataset containing the values of the distribution at each
        point. Each value must be associated with a corresponding point in the
        `position` property.
        """
        return self._data

    @property
    def resolution(self):
        """float: Resolution of the dataset, in units of points / mm."""
        return self._resolution

    @property
    def translation(self):
        """ndarray: Represents a shift applied to the distribution from its
        original position. For example, a shift of 1 mm in the x-direction and
        -2 mm in the y-direction would be np.array([1, -2]).
        """
        return self._translation

    @translation.setter
    def translation(self, translation):
        """ndarray: Represents a shift applied to the distribution from its
        original position. For example, a shift of 1 mm in the x-direction and
        -2 mm in the y-direction would be np.array([1, -2]).
        """
        self._translation = translation

    @property
    def rotation(self):
        """float: Represents a rotation of the dataset from its original
        location, in units of radians.
        """
        return self._rotation

    @rotation.setter
    def rotation(self, rad):
        self._rotation = rad

    @property
    def position(self):
        """ndarray: The spatial position of each point in the distribution.
        For a 2-dimensional distribution of size M x N, the position has
        shape (2, M, N), with the first dimension corresponding to the x and y
        coordinates of each point.
        """
        # Rotate
        # TODO: vectorise this algorithm
        pos = copy.deepcopy(self._position)
        if self._rotation:
            rot_matrix = (
                np.array([[np.cos(self._rotation), -np.sin(self._rotation)],
                          [np.sin(self._rotation), np.cos(self._rotation)]])
            )
            for i in range(self._position.shape[1]):
                for j in range(self._position.shape[2]):
                    pos[:, i, j] = np.matmul(rot_matrix, pos[:, i, j])
        # Translate
        pos = pos + self._translation[:, np.newaxis, np.newaxis]
        return pos

    @position.setter
    def position(self, pos):
        self._position = pos

    def scale_grid(self, reference_distribution=None, new_resolution=None):
        """Increase or decrease the resolution of a distribution.

        Scales the resolution of the distribution by first interpolating the
        data values and then resampling it at the positions specified in the
        reference distribution. If the `scale` parameter is specified, the
        number of points in each dimension will be increased `scale` times.
        Currently, only regularly spaced gridded data is supported, due to the
        use of RecBivarateSpline for interpolation. See the
        scipy.interpolate.RectBivariateSpline documentation for more
        information.

        Parameters
        ----------
        reference_distribution : Distribution
            Interpolate new data values at the positions in this reference
            distribution.
        new_resolution : float
            The desired resolution of the scaled distribution, in pixels / mm.
            If None, the scaled distribution will have the same resolution as
            reference_distribution.

        Returns
        -------
        Distribution
            A new distrubution with points at locations given by the reference
            distribution, scaled appropriately.
        """
        if reference_distribution is None:
            reference_distribution = self

        if new_resolution is None:
            new_resolution = reference_distribution.resolution
            scale = 1
        else:
            scale = new_resolution / reference_distribution.resolution

        # Generate new positions by scaling the positions given in
        # reference_distribution
        x = reference_distribution.position[0, 0, :]
        y = reference_distribution.position[1, :, 0]
        xscale = scale * (len(x) - 1) + 1
        yscale = scale * (len(y) - 1) + 1
        x_new = np.linspace(x[0], x[-1], xscale)
        y_new = np.linspace(y[0], y[-1], yscale)

        # Set up linear interpolation function for original grid
        interp_spline = interpolate.RectBivariateSpline(
            self.position[0, 0, :],
            self.position[1, :, 0],
            self.data.T,
            kx=1,
            ky=1
        )

        # Find data points at new positions
        new_data = interp_spline(x_new, y_new)
        xx, yy = np.meshgrid(x_new, y_new)
        new_distribution = Distribution(
            new_data.T,
            resolution=new_resolution,
            position=np.stack([xx, yy])
        )
        return new_distribution
