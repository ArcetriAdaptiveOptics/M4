#!/usr/bin/env python
import unittest
from m4.devices.parabola_slider import OpcUaParabolaSlider
from unittest.mock import MagicMock, Mock
from m4.configuration.ott_parameters import OpcUaParameters


class TestOpcUaParabolaSlider(unittest.TestCase):

    def setUp(self):
        self.opc = MagicMock()
        self.slider = OpcUaParabolaSlider(self.opc)

    def testCallsOpcOnGetPosition(self):
        self.opc.get_position = Mock(return_value=123.4)
        pos = self.slider.getPosition()
        self.opc.get_position.assert_called_once_with(OpcUaParameters.ST)
        self.assertAlmostEqual(pos, 123.4)

    def testCallsOpcOnSetPosition(self):
        self.opc.get_position = Mock(return_value=3.14)
        pos = self.slider.setPosition(3.14)
        self.opc.set_target_position.assert_called_once_with(
            OpcUaParameters.ST, 3.14)
        self.opc.move_object.assert_called_once_with(OpcUaParameters.ST)
        self.opc.wait_for_stop.assert_called_once_with(OpcUaParameters.ST)
        self.assertAlmostEqual(pos, 3.14)

    def testMoveRelativeToCurrentPosition(self):

        self._mocked_pos = {}

        def side_effect_get(key):
            return self._mocked_pos[key]

        def side_effect_set(key, val):
            self._mocked_pos[key] = val
            print("%s" % self._mocked_pos)

        self.opc.get_position = Mock(side_effect=side_effect_get)
        self.opc.set_target_position = Mock(side_effect=side_effect_set)

        init_pos = 42
        self.slider.setPosition(init_pos)
        curr_pos = self.slider.getPosition()

        off = 19.3
        fin_pos = self.slider.setPosition(curr_pos + off)

        self.assertAlmostEqual(init_pos + off, fin_pos)
