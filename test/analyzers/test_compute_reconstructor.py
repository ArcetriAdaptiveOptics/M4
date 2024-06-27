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
        random_matrix = [[random.random() for i in range(100)]
                         for j in range(100)]
        self._random_matrix_dp = np.linalg.inv(
            np.array(random_matrix).T @ np.array(random_matrix)
        )

    @mock.patch("m4.analyzers.analyzer_iffunctions.AnalyzerIFF")
    @mock.patch("matplotlib.pyplot")
    @mock.patch("matplotlib.pyplot.show")
    def test_compute_reconstructor(self, m_analyzer, m_plot, m_show):
        m_analyzer.getInteractionMatrix.return_value = self._random_matrix_dp
        m_plot.return_value = {'y': 1, 'x': 1}
        rec = ComputeReconstructor(m_analyzer)
        _ = rec.run(Interactive=False)
        _ = rec.run(Interactive=False, sv_threshold=3.9e-2)
        np.testing.assert_array_almost_equal_nulp(
            rec._filtered_sv[-7:], np.zeros(7))

    def tearDown(self):
        self._random_matrix_dp = None
        pass
