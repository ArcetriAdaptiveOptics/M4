'''
Authors
  - C. Selmi: written in 2020
'''
import os
import unittest
from types import SimpleNamespace
from unittest.mock import patch, MagicMock
from test.helper_test_library import testDataRootDir
from m4.devices.accelerometers import ZmqAccelerometers

TT = '20210519_160814'


class TestZmqAccelerometers(unittest.TestCase):

    def setUp(self):
        self.acc = ZmqAccelerometers()
        self.data_folder = os.path.join(
            testDataRootDir(), 'base', 'M4Data',
            'OPTData', 'AccelerometersData')

    @patch('m4.devices.accelerometers._time', MagicMock())
    @patch('m4.devices.accelerometers._ntn', return_value=TT)
    @patch('m4.devices.accelerometers._zmq')
    def testForDataAcquisition(self, zmq_mock, ntn_mock):
        # the real zmq socket is replaced: the server "replies" with one byte
        socket = zmq_mock.Context.return_value.socket.return_value
        socket.recv.return_value = b'0'
        fake_params = SimpleNamespace(
            accelerometers_server='tcp://fake:1234',
            accelerometers_data_folder=self.data_folder)
        fake_folders = SimpleNamespace(
            ACCELEROMETERS_ROOT_FOLDER='/fake/root')
        self.acc._acc = MagicMock()

        with patch('m4.devices.accelerometers._opcua', fake_params), \
                patch('m4.devices.accelerometers._fn', fake_folders):
            tt = self.acc.acquireData(3)

        self.assertEqual(tt, TT)
        socket.connect.assert_called_once_with('tcp://fake:1234')
        socket.send_string.assert_called_once_with('START 3')
        socket.disconnect.assert_called_once_with('tcp://fake:1234')
        # the last file of the data folder is converted to the final place
        self.acc._acc.convertAndSaveData.assert_called_once_with(
            os.path.join(self.data_folder, TT + '.h5'),
            os.path.join('/fake/root', TT + '.h5'))

    @patch('m4.devices.accelerometers._time', MagicMock())
    @patch('m4.devices.accelerometers._zmq')
    def testNoReplyFromSocketRaises(self, zmq_mock):
        socket = zmq_mock.Context.return_value.socket.return_value
        socket.recv.side_effect = Exception('timeout')
        with self.assertRaises(OSError):
            self.acc.acquireData(1)
