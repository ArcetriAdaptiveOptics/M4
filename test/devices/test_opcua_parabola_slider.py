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
        self.opc.get_positions = Mock(return_value=123.4)
        pos = self.slider.getPosition()
        self.opc.get_positions.assert_called_once_with(OpcUaParameters.ST)
        self.assertAlmostEqual(pos, 123.4)

    def testCallsOpcOnSetPosition(self):
        self.opc.get_positions = Mock(return_value=3.14)
        # setPosition prints the reached position and returns nothing
        self.slider.setPosition(3.14)
        self.opc.set_position.assert_called_once_with(OpcUaParameters.ST, 3.14)
        self.opc.wait_for_stop.assert_called_once_with(OpcUaParameters.ST)
        self.opc.get_positions.assert_called_once_with(OpcUaParameters.ST)

    def testMoveRelativeToCurrentPosition(self):
        self._mocked_pos = {}

        def side_effect_get(key):
            return self._mocked_pos[key]

        def side_effect_set(key, val):
            self._mocked_pos[key] = val

        self.opc.get_positions = Mock(side_effect=side_effect_get)
        self.opc.set_position = Mock(side_effect=side_effect_set)

        init_pos = 42
        self.slider.setPosition(init_pos)
        curr_pos = self.slider.getPosition()

        off = 19.3
        self.slider.setPosition(curr_pos + off)
        fin_pos = self.slider.getPosition()

        self.assertAlmostEqual(init_pos + off, fin_pos)

    def testGetPositionInMeters(self):
        self.opc.get_positions = Mock(return_value=1500.0)
        self.assertAlmostEqual(self.slider.getPositionInM(), 1.5)

    def testOutOfRangeSlideRaises(self):
        with self.assertRaises(ValueError):
            self.slider.setPosition(OpcUaParameters.max_slide + 1)
        with self.assertRaises(ValueError):
            self.slider.setPosition(OpcUaParameters.min_slide - 1)
        self.opc.set_position.assert_not_called()
