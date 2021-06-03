'''
Authors
  - C. Selmi: written in 2020
'''
import os
import unittest
from test.test_helper import testDataRootDir
from m4.devices.accelerometers import ZmqAccelerometes
import mock


class TestZmqAccelerometers(unittest.TestCase):

    def setUp(self):
        self.acc = ZmqAccelerometes()

    @mock.patch('m4.devices.accelerometers.fold_name', autospect=True)
    @mock.patch('m4.devices.accelerometers.OpcUaParameters', autospect=True)
    @mock.patch('zmq.Context', autospect=True)
    def testForDataAcquisition(self, mock_final_fold_name, mock_start_fold_name, zmq_mock):
        want_acc_root_folder = os.path.join(
            testDataRootDir(), 'base', 'M4Data',
            'OPTData', 'AccelerometersData')
        mock_final_fold_name.ACC_ROOT_FOLDER = want_acc_root_folder
        mock_start_fold_name.accelerometers_data_folder = want_acc_root_folder
        print(mock_start_fold_name.accelerometers_data_folder)

        acc = ZmqAccelerometes()
        #name = acc.acquireData()
