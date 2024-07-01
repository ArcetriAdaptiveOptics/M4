# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 11:54:14 2024

Author(s): 
    - P. Ferraiuolo
"""
# Test iff_acquisition_preparation
import unittest
# from unittest.mock import Mock
import numpy as np
from m4.dmutils import iff_acquisition_preparation as ifa

class TestIFFCapturePreparation(unittest.TestCase):

    def setUp(self):
        # Fake dm for initialization
        self.mock_dm = type('MockDM', (object,), {
            "mirrorModes": np.array([[1, 2], [3, 4]]),
            "nActs": 2
        })()
        self.iff_capture_preparation = ifa.IFFCapturePreparation(self.mock_dm)

    def test_initialization(self):
        self.assertEqual(self.iff_capture_preparation._NActs, 2)
        np.testing.assert_array_equal(self.iff_capture_preparation._mirrorModes, np.array([[1, 2], [3, 4]]))
        self.assertIsNone(self.iff_capture_preparation._cmdMatrix)
        self.assertEqual(self.iff_capture_preparation._modalBase, self.iff_capture_preparation._mirrorModes)
        
    def test_createTimedCmdHistory(self):
        modes_list = np.array([0,1,2])
        amplitude = 1e-6
        template = [1,-1,1]
        self.iff_capture_preparation.createTimedCmdHistory(modesList=modes_list, 
                                                           modesAmp=amplitude, 
                                                           template=template)
        np.testing.assert_
        # Verificare il risultato atteso

        self.assertIsNotNone(self.iff_capture_preparation.cmdMatHistory)
