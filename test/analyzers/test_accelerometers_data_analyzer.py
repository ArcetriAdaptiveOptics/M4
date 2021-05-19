import unittest
import os
from m4.analyzers.accelerometers_data_analyzer import AccelerometersDataAnalyzer
import mock
from test.test_helper import testDataRootDir


class TestAccelerometersDataAnalyzer(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @mock.patch('m4.analyzers.accelerometers_data_analyzer.fold_name', autospect=True)
    def testName(self, mock_fold_name):
        want_acc_root_folder = os.path.join(
            testDataRootDir(), 'base', 'M4Data',
            'OPTData', 'AccelerometersData')
        mock_fold_name.ACC_ROOT_FOLDER = want_acc_root_folder

        tt = '20210519_160814'
        self.ana = AccelerometersDataAnalyzer(tt)
        pass

