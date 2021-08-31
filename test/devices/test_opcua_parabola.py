#!/usr/bin/env python
import unittest
import numpy as np
from m4.devices.parabola import OpcUaParabola
import mock
from m4.configuration.ott_parameters import OpcUaParameters
from mock.mock import call
from numpy import testing


class TestOpcUaParabola(unittest.TestCase):

    def setUp(self):
        self._opc = mock.MagicMock()
        self._parabola = OpcUaParabola(self._opc)

    def testGetPositionReturnSixTuple(self):
        pos = self._parabola.getPosition()
        self.assertTrue(pos.shape, (6,))

    def testCallOpcUaOnGetPosition(self):
        piston = -11
        tip = 1e-5
        tilt = 0.003

        values = {OpcUaParameters.PAR_PISTON: piston,
                  OpcUaParameters.PAR_TILT: tilt,
                  OpcUaParameters.PAR_TIP: tip}

        def side_effect(arg):
            return values[arg]

        self._opc.get_position = mock.Mock(side_effect=side_effect)

        pos = self._parabola.getPosition()
        calls = [call(OpcUaParameters.PAR_PISTON),
                 call(OpcUaParameters.PAR_TILT),
                 call(OpcUaParameters.PAR_TIP)]
        self._opc.get_position.assert_has_calls(calls, any_order=True)
        testing.assert_almost_equal(
            pos,
            np.array([0, 0, piston, tip, tilt, 0]))

    def testCallOpcUaOnSetPosition(self):
        piston = -123
        tip = 3.13
        tilt = -12e-7
        pos = np.array([0, 0, piston, tip, tilt, 0])

        self._values = {}

        def side_effect_set(dof, value):
            self._values[dof] = value

        def side_effect_get(arg):
            return self._values[arg]

        self._opc.get_position = mock.Mock(side_effect=side_effect_get)
        self._opc.set_target_position = mock.Mock(side_effect=side_effect_set)

        newpos = self._parabola.setPosition(pos)

        calls = [call(OpcUaParameters.PAR_PISTON, piston),
                 call(OpcUaParameters.PAR_TILT, tilt),
                 call(OpcUaParameters.PAR_TIP, tip)]
        self._opc.set_target_position.assert_has_calls(calls, any_order=True)
        self._opc.move_object.assert_called_once_with(OpcUaParameters.PAR_KIN)
        self._opc.wait_for_stop.assert_called_once_with(OpcUaParameters.PAR_KIN)

        testing.assert_almost_equal(pos, newpos)

