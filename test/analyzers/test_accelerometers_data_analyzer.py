'''
Authors
  - C. Selmi:  written in 2021
'''
import unittest
import os
from m4.analyzers.accelerometers_data_analyzer import AccelerometersDataAnalyzer
from unittest.mock import patch
from test.helper_test_library import testDataRootDir


class TestAccelerometersDataAnalyzer(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @patch('m4.type.accelerometers_data.fold_name', unsafe=True)
    def testReadingAndAnalysis(self, mock_fold_name):
        want_acc_root_folder = os.path.join(
            testDataRootDir(), 'base', 'M4Data',
            'OPTData', 'AccelerometersData')
        mock_fold_name.ACC_ROOT_FOLDER = want_acc_root_folder

        tt = '20210519_160814'
        self.ana = AccelerometersDataAnalyzer(tt)
        #self.ana.readAndShow()
