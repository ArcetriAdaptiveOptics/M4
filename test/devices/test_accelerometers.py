'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
from m4.devices.accelerometers import ZmqAccelerometes
import mock
from m4.configuration.ott_parameters import OpcUaParameters


class TestZmqAccelerometers(unittest.TestCase):

    def setUp(self):
        self.acc = ZmqAccelerometes()


    @mock.patch('zmq', autospect=True)
    def testForDataAcquisition(self, zmq_mock):
        pass
