'''
Authors
  - C. Selmi: written in 2021
'''
from numpy import testing
import numpy as np
import unittest
from unittest.mock import MagicMock, Mock
from m4.devices.temperature_sensors import OpcUaTemperatureSensors
from m4.configuration.ott_parameters import OpcUaParameters

class TestOpcTemperatureSensor(unittest.TestCase):

    def setUp(self):
        self.opc = MagicMock()
        self.temp = OpcUaTemperatureSensors(self.opc)

    def testCallsOpcOnGetTemperature(self):
        self.opc.get_temperature_vector = Mock(return_value=np.ones(24))
        ret = self.temp.getTemperature()
        testing.assert_almost_equal(ret,
                                    np.ones(OpcUaParameters.num_PT_sensor))

    def testGetPositionReturnArray(self):
        self.opc.get_temperature_vector = Mock(return_value=np.ones(24))
        ret = self.temp.getTemperature()
        message = "Get position doesn't return a vector of shape 24."
        self.assertEqual(ret.shape, (24,), message)
