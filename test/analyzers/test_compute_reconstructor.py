'''
Authors
  - C. Selmi:  written in 2022
'''
import mock
import numpy as np
import matplotlib.pyplot as plt
from m4.analyzers.analyzer_iffunctions import AnalyzerIFF
from m4.analyzers.compute_reconstructor import ComputeReconstructor
import random
import unittest


class TestComputeReconstructor(unittest.TestCase):

    def setUp(self):
        random.seed(0)
        random_matrix = [
            [random.random() for i in range(100)] for j in range(1000)]
        self._random_matrix_dp = np.linalg.inv(
            np.array(random_matrix).T @ np.array(random_matrix))

    @mock.patch('m4.analyzers.analyzer_iffunctions.AnalyzerIFF')
    def test_compute_reconstructor(self, mock_analyzer):
        mock_analyzer.getInteractionMatrix.return_value = self._random_matrix_dp
        rec = ComputeReconstructor(mock_analyzer)
        _ = rec.run(Interactive=False)
        _ = rec.run(Interactive=False, sv_threshold=8e-3)
        print(rec._intMat_S)
        print(rec._threshold)
        print(rec._filtered_sv)
        self.assertAlmostEqual(rec._intMat_S, rec._intMat_S*0)

    def tearDown(self):
        pass
