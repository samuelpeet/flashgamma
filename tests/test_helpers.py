import os
import pytest
import numpy as np
from flashgamma import load_file, create_distance_kernel


class TestLoadFile:
    def test_load_file_arccheck(self):
        pass

    def test_load_file_snc(self):
        pass

    def test_load_file_dicom(self, tmpdir):
        f = tmpdir.join("blank_test_file.dcm")
        f.write("foo")
        with pytest.raises(NotImplementedError):
            load_file(f, "dicom")

    def test_load_file_tiff(self, tmpdir):
        f = tmpdir.join("blank_test_file.tiff")
        f.write("foo")
        with pytest.raises(NotImplementedError):
            load_file(f, "tiff")

    def test_load_file_unknown(self, tmpdir):
        f = tmpdir.join("blank_test_file.png")
        f.write("foo")
        with pytest.raises(NotImplementedError):
            load_file(f, "png")

class TestCreateDistanceKernal:
    def test_create_distance_kernel(self):
        dta = 2
        resolution = 0.5
        kernel = create_distance_kernel(dta, resolution)
        correct = np.array([[2.5e8, 1.0e0, 2.5e8],
                            [1.0e0, 0.0e0, 1.0e0],
                            [2.5e8, 1.0e0, 2.5e8]])
        np.testing.assert_equal(kernel, correct)
