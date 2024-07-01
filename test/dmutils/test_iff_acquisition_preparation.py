'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import numpy as np
import os
import mock
from test.helper_test_library import testDataRootDir
from m4.utils import optical_calibration
from dmutils.iff_acquisition_preparation import IffAcquisitionPreparation

class TestIffAcquisitionPreparation(unittest.TestCase):


    def setUp(self):


    def tearDown(self):
        del(self.cal, self.ott)


