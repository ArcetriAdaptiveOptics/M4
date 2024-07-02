# -*- coding: utf-8 -*-
"""
Created on July 2024

Author(s): 
    - P. Ferraiuolo
"""
# Test iff_acquisition_preparation
import unittest
from unittest import mock
from parameterized import parameterized
import numpy as np
from m4.dmutils import iff_acquisition_preparation as ifa

class TestIFFCapturePreparation(unittest.TestCase):

    def setUp(self):
        # Fake dm for initialization
        mockMirrorModes = np.zeros((10,10), int)
        np.fill_diagonal(mockMirrorModes, 5)
        self.mock_dm = type('MockDM', (object,), {
            "mirrorModes": mockMirrorModes,
            "nActs": 10
        })()
        self.iff_capture_preparation = ifa.IFFCapturePreparation(self.mock_dm)

    def test_initialization(self):
        self.assertEqual(self.iff_capture_preparation._NActs, 10)
        np.testing.assert_array_equal(self.iff_capture_preparation.mirrorModes, 
                                      np.array([[5, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                             [0, 5, 0, 0, 0, 0, 0, 0, 0, 0],
                                             [0, 0, 5, 0, 0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 5, 0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 5, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 5, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0, 5, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0, 0, 5, 0, 0],
                                             [0, 0, 0, 0, 0, 0, 0, 0, 5, 0],
                                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 5]]))
        self.assertIsNone(self.iff_capture_preparation._cmdMatrix)
        np.testing.assert_array_equal(self.iff_capture_preparation._modalBase, 
                                      self.iff_capture_preparation.mirrorModes)
        
    @parameterized.expand(['mirror','zonal','hadamard'])
    def test_updateModalBase(self, base_name):
        self.iff_capture_preparation._updateModalBase(mbasename = base_name)
        self.assertIsNotNone(self.iff_capture_preparation._modalBase)
        self.assertEqual(self.iff_capture_preparation.modalBaseId, base_name)
        self.assertEqual(self.iff_capture_preparation._modalBase.shape, 
                         (self.iff_capture_preparation._NActs, self.iff_capture_preparation._NActs))
        
    @mock.patch('m4.configuration.read_iffconfig', new='mock_read_iffconfig')
    def test_createTriggerPadding(self, mock_read_iffconfig):
        mock_read_iffconfig.getConfig('TRIGGER').return_value = {
            'modaleBase': 'mirror', 
            'zeros': 5,
            'modes': 5,
            'amplitude': 1e-6,
            'template': [1]
            }
        self.iff_capture_preparation._createTriggerPadding()
        self.assertIsNotNone(self.iff_capture_preparation.triggPadCmdHist)
        self.assertEqual(self.iff_capture_preparation.triggPadCmdHist.shape, 
                         (self.iff_capture_preparation._NActs, 6))
        self.assertIsNotNone(self.iff_capture_preparation.modalBaseId)
        self.assertEqual(self.iff_capture_preparation.modalBaseId, 'mirror')
        
    @mock.patch('m4.configuration.read_iffconfig', new='mock_read_iffconfig')
    def test_createRegistrationPattern(self, mock_read_iffconfig):
        mock_read_iffconfig.getConfig('REGISTRATION').return_value = {
            'modaleBase': 'zonal',
            'zeros': 1,
            'modes': [110,350,600],
            'amplitude': 1e-6,
            'template': [1,-1]
            }
        self.iff_capture_preparation._createRegistrationPattern()
        self.assertIsNotNone(self.iff_capture_preparation.regPadCmdHist)
        self.assertEqual(self.iff_capture_preparation.regPadCmdHist.shape, 
                         (self.iff_capture_preparation._NActs, 7))
        self.assertIsNotNone(self.iff_capture_preparation.modalBaseId)
        self.assertEqual(self.iff_capture_preparation.modalBaseId, 'mirror')
                
    def test_createTimedCmdHistory(self, mock_createAuxCmdHistory, mock_createCmdMatrixHistory, mock_read_iffconfig):
        mock_createAuxCmdHistory.return_value = np.fill((self.iff_capture_preparation._NActs, 13),
                                                        1)
        mock_createCmdMatrixHistory.return_value = np.fill((self.iff_capture_preparation._NActs, 9),
                                                           2)
        mock_read_iffconfig.getTiming.return_value = 1
        # modes_list = np.array([0,1,2])
        # amplitude = 2
        # template = [1,-1,1]
        self.iff_capture_preparation.createTimedCmdHistory()
        self.assertIsNotNone(self.iff_capture_preparation.cmdMatHistory)
        self.assertEqual(self.iff_capture_preparation.timedCmdHistory.shape, 
                         (self.iff_capture_preparation._NActs, 22))
