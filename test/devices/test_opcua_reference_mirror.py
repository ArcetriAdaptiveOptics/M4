'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
import numpy as np
from m4.devices.reference_mirror import OpcUaReferenceMirror
from unittest.mock import Mock, MagicMock
from m4.configuration.ott_parameters import OpcUaParameters
from numpy import testing


class TestOpcUaReferenceMirror(unittest.TestCase):

    def setUp(self):
        self._opc = MagicMock()
        self._rm = OpcUaReferenceMirror(self._opc)
        # the device is a tripod: the controller exposes [piston, tip, tilt]
        self._opc.get_act_positions = Mock(return_value=[0.0, 0.0, 0.0])

    def testGetPositionReturnSixTuple(self):
        pos = self._rm.getPosition()
        self.assertEqual(pos.shape, (6,),
                         "Get position doesn't return a vector of shape 6.")

    def testCallOpcUaOnGetPosition(self):
        piston = 1e-16
        tip = 1e-6
        tilt = 0.007
        self._opc.get_act_positions = Mock(return_value=[piston, tip, tilt])

        pos = self._rm.getPosition()

        self._opc.get_act_positions.assert_called_once_with(OpcUaParameters.RM)
        testing.assert_almost_equal(
            pos, np.array([0, 0, piston, tip, tilt, 0]))

    def testCallOpcUaOnSetPosition(self):
        piston = -7
        tip = 30.23
        tilt = -12e-3
        pos = np.array([0, 0, piston, tip, tilt, 0])

        self._values = {}

        def side_effect_set(dev, values):
            self._values[dev] = list(values)

        def side_effect_get(dev):
            return self._values[dev]

        self._opc.get_act_positions = Mock(side_effect=side_effect_get)
        self._opc.set_act_positions = Mock(side_effect=side_effect_set)

        # setPosition prints the reached position and returns nothing
        self._rm.setPosition(pos)

        self._opc.set_act_positions.assert_called_once()
        args = self._opc.set_act_positions.call_args.args
        self.assertEqual(args[0], OpcUaParameters.RM)
        testing.assert_almost_equal(args[1], [piston, tip, tilt])
        self._opc.wait_for_stop.assert_called_once_with(OpcUaParameters.RM)

        testing.assert_almost_equal(pos, self._rm.getPosition())
