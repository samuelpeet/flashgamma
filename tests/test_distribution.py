import pytest
import numpy as np
from flashgamma import Distribution

@pytest.fixture()
def simple_dist():
    data = np.zeros((11, 11))
    data[3:8, 3:8] = 2.0
    return Distribution(data)   

@pytest.fixture()
def simple_dist_with_position():
    data = np.zeros((11, 11))
    data[3:8, 3:8] = 2.0
    x = np.linspace(-4, 6, 11)
    y = np.linspace(-3, 7, 11)
    xx, yy = np.meshgrid(x, y)
    position = np.stack([xx, yy])
    return Distribution(data, position=position)

class TestConstructor:
    def test_constructor(self, simple_dist):
        x = np.linspace(-5, 5, 11)
        y = np.linspace(-5, 5, 11)
        xx, yy = np.meshgrid(x, y)
        correct_pos = np.stack([xx, yy])
        np.testing.assert_array_equal(simple_dist.position, correct_pos)

        correct_data = np.zeros((11, 11))
        correct_data[3:8, 3:8] = 2.0
        np.testing.assert_array_equal(simple_dist.data, correct_data)
    
    def test_constructor_with_position(self, simple_dist_with_position):
        x = np.linspace(-4, 6, 11)
        y = np.linspace(-3, 7, 11)
        xx, yy = np.meshgrid(x, y)
        correct = np.stack([xx, yy])
        np.testing.assert_array_equal(
            simple_dist_with_position.position, correct
        )

class TestTranslation:
    def test_translation(self, simple_dist):
        x = np.linspace(-4, 6, 11)
        y = np.linspace(-3, 7, 11)
        xx, yy = np.meshgrid(x, y)
        correct = np.stack([xx, yy])
        simple_dist.translation = np.array([1, 2])
        np.testing.assert_array_equal(
            simple_dist.translation,
            np.array([1, 2])
        )
        np.testing.assert_array_equal(
            simple_dist.position, correct
        )

class TestRotation:
    def test_rotation(self, simple_dist):
        x = np.linspace(-5, 5, 11)
        y = np.linspace(-5, 5, 11)
        xx, yy = np.meshgrid(x, y)
        pos_rotated = np.stack([xx, yy])       
        pos_rotated = np.array([
            np.rot90(pos_rotated[0]),
            np.rot90(pos_rotated[1])
        ])

        simple_dist.rotation = np.pi / 2
        assert simple_dist.rotation == np.pi / 2   
        np.testing.assert_array_almost_equal(
            simple_dist.position, pos_rotated
        )

class TestScaleGrid:
    def test_scale_grid(self, simple_dist):
        # Recreate simple_dist with double resolution
        data = np.zeros((21, 21))
        data[6:16, 6:16] = 2.0
        x = np.linspace(-5, 5, 21)
        y = np.linspace(-5, 5, 21)
        xx, yy = np.meshgrid(x, y)
        position = np.stack([xx, yy])
        simple_dist_double_res = Distribution(data, position=position)

        simple_dist = simple_dist.scale_grid(new_resolution=2)
        np.testing.assert_array_almost_equal(
            simple_dist.position,
            simple_dist_double_res.position
        )
