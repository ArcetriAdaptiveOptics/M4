#!/usr/bin/env python
import unittest
import numpy as np
from m4.devices.parabola import OpcUaParabola
from unittest.mock import Mock, MagicMock
from m4.configuration.ott_parameters import OpcUaParameters
from numpy import testing


class TestOpcUaParabola(unittest.TestCase):

    def setUp(self):
        self._opc = MagicMock()
        self._parabola = OpcUaParabola(self._opc)
        # the device is a tripod: the controller exposes [piston, tip, tilt]
        self._opc.get_act_positions = Mock(return_value=[0.0, 0.0, 0.0])

    def testGetPositionReturnSixTuple(self):
        pos = self._parabola.getPosition()
        self.assertEqual(pos.shape, (6,),
                         "Get position doesn't return a vector of shape 6.")

    def testCallOpcUaOnGetPosition(self):
        piston = -11
        tip = 1e-5
        tilt = 0.003
        self._opc.get_act_positions = Mock(return_value=[piston, tip, tilt])

        pos = self._parabola.getPosition()

        self._opc.get_act_positions.assert_called_once_with(OpcUaParameters.PAR)
        testing.assert_almost_equal(
            pos, np.array([0, 0, piston, tip, tilt, 0]))

    def testCallOpcUaOnSetPosition(self):
        piston = -123
        tip = 3.13
        tilt = -12e-7
        pos = np.array([0, 0, piston, tip, tilt, 0])

        self._values = {}

        def side_effect_set(dev, values):
            self._values[dev] = list(values)

        def side_effect_get(dev):
            return self._values[dev]

        self._opc.get_act_positions = Mock(side_effect=side_effect_get)
        self._opc.set_act_positions = Mock(side_effect=side_effect_set)

        # setPosition prints the reached position and returns nothing
        self._parabola.setPosition(pos)

        self._opc.set_act_positions.assert_called_once()
        args = self._opc.set_act_positions.call_args.args
        self.assertEqual(args[0], OpcUaParameters.PAR)
        testing.assert_almost_equal(args[1], [piston, tip, tilt])
        self._opc.wait_for_stop.assert_called_once_with(OpcUaParameters.PAR)

        testing.assert_almost_equal(pos, self._parabola.getPosition())
