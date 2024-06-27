"""
Authors
  - C. Selmi:  written in 2022
"""

import mock
import numpy as np
from m4.analyzers.compute_reconstructor import ComputeReconstructor
import random
import unittest


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
        _ = self._cr.run(Interactive=False, sv_threshold=11)
        np.testing.assert_array_almost_equal_nulp(
            self._cr._filtered_sv[-11:], np.zeros(11))
        # _ = self._cr.run(Interactive=True)

    @mock.patch("astropy.io.fits")
    def test_load_reconstructor(self, m_fits):
        m_fits.open.return_value = self._intMatCube
        tn = "20220101_000000"
        self._cr = ComputeReconstructor.loadIntMatFromFolder(tn)

    def tearDown(self):
        self._random_matrix_dp = None
        pass
