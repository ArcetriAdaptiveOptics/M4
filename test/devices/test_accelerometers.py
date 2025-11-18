'''
Authors
  - C. Selmi: written in 2020
'''
import os
import unittest
from test.helper_test_library import testDataRootDir
from m4.devices.accelerometers import ZmqAccelerometers
from unittest.mock import patch
from m4.type.accelerometers_data import AccelerometersData


class TestZmqAccelerometers(unittest.TestCase):

    def setUp(self):
        self.acc = ZmqAccelerometers()

    @patch('m4.devices.accelerometers.fold_name', unsafe=True)
    @patch('m4.devices.accelerometers.OpcUaParameters', unsafe=True)
    @patch('zmq.Context', unsafe=True)
    def testForDataAcquisition(self, mock_final_fold_name, mock_start_fold_name, zmq_mock):
        want_acc_root_folder = os.path.join(
            testDataRootDir(), 'base', 'M4Data',
            'OPTData', 'AccelerometersData')
        mock_final_fold_name.ACC_ROOT_FOLDER = want_acc_root_folder
        mock_start_fold_name.accelerometers_data_folder = want_acc_root_folder
        print(mock_start_fold_name.accelerometers_data_folder)

        acc = ZmqAccelerometers()
        #name = acc.acquireData()

        tt = '20210519_160814.h5'
        acc = AccelerometersData()
        start = os.path.join(mock_start_fold_name.accelerometers_data_folder, tt)
        final_destination = os.path.join(mock_final_fold_name.ACC_ROOT_FOLDER, tt)
        #acc.convertAndSaveData(start, final_destination)

