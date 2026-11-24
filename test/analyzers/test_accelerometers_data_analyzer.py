'''
Authors
  - C. Selmi:  written in 2021
'''
import unittest
import os
from types import SimpleNamespace
from m4.analyzers.accelerometers_data_analyzer import AccelerometersDataAnalyzer
from unittest.mock import patch
from test.helper_test_library import testDataRootDir


class TestAccelerometersDataAnalyzer(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testReadingAndAnalysis(self):
        want_acc_root_folder = os.path.join(
            testDataRootDir(), 'base', 'M4Data',
            'OPTData', 'AccelerometersData')
        fake_folders = SimpleNamespace(
            ACCELEROMETERS_ROOT_FOLDER=want_acc_root_folder)

        tt = '20210519_160814'
        with patch('m4.type.accelerometers_data.fold_name', fake_folders):
            self.ana = AccelerometersDataAnalyzer(tt)
        self.assertEqual(self.ana.tt, tt)
        self.assertEqual(
            self.ana._h5_file_path,
            os.path.join(want_acc_root_folder, tt + '.h5'))
        #self.ana.readAndShow()
