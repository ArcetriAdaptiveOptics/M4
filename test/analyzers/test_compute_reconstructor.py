"""
Authors
  - C. Selmi:  written in 2022
"""

import mock
import numpy as np
from m4.analyzers.compute_reconstructor import ComputeReconstructor
import random
import unittest


class fake_hdu():

    def __init__(self, data):
        self.data = data

    def __getitem__(self, key):
        return self

    def close(self):
        pass


class TestComputeReconstructor(unittest.TestCase):

    def setUp(self):
        random.seed(0)
        self._intMatCube = []
        for k in range(50):
            random_matrix = [[random.random() for i in range(100)]
                             for j in range(100)]
            x, y = np.meshgrid(np.arange(100)-50+3*np.random.random(),
                               np.arange(100)-50+3*np.random.random())
            mask = (x**2 + y**2) > 25**2
            self._intMatCube.append(
                np.ma.masked_array(random_matrix, mask=mask))
        # self._intMatCube = np.array(self._intMatCube)
        self._intMatCube = np.ma.dstack(self._intMatCube)

    def test_compute_reconstructor(self):
        self._cr = ComputeReconstructor(self._intMatCube)
        _ = self._cr.run(Interactive=False)

        self._cr = ComputeReconstructor(self._intMatCube)
        _ = self._cr.run(Interactive=False, sv_threshold=11)

        print(self._cr._filtered_sv[-10:])
        np.testing.assert_array_almost_equal_nulp(
            self._cr._filtered_sv[-10:], np.zeros(10))
        # _ = self._cr.run(Interactive=True)

    @mock.patch("os.path.join")
    @mock.patch("astropy.io.fits.open")
    def test_load_reconstructor(self, m_fits, m_join):
        m_join.return_value = "20220101_000000"
        m_fits.side_effect =\
            [fake_hdu(self._intMatCube.copy()), fake_hdu(
                self._intMatCube.copy())]  # finire questo test!
        tn = "20220101_000000"
        self._cr = ComputeReconstructor.loadIntMatFromFolder(tn)
        _ = self._cr.run(Interactive=False, sv_threshold=11)

        print(self._cr._filtered_sv[-10:])
        np.testing.assert_array_almost_equal_nulp(
            self._cr._filtered_sv[-10:], np.zeros(10))

    def tearDown(self):
        self._random_matrix_dp = None
        self._cr = None
        pass
