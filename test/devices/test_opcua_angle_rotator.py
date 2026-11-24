'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
from unittest.mock import MagicMock, Mock
from m4.devices.angle_rotator import OpcUaAngleRotator
from m4.configuration.ott_parameters import OpcUaParameters


class TestOpcUaAngleRotator(unittest.TestCase):

    def setUp(self):
        self.opc = MagicMock()
        self.angle = OpcUaAngleRotator(self.opc)

    def testCallsOpcOnGetPosition(self):
        self.opc.get_positions = Mock(return_value=60)
        pos = self.angle.getPosition()
        self.opc.get_positions.assert_called_once_with(OpcUaParameters.RA)
        self.assertAlmostEqual(pos, 60)

    def testCallsOpcOnSetPosition(self):
        self.opc.get_positions = Mock(return_value=3.14)
        # setPosition prints the reached position and returns nothing
        self.angle.setPosition(3.14)
        self.opc.set_position.assert_called_once_with(OpcUaParameters.RA, 3.14)
        self.opc.wait_for_stop.assert_called_once_with(OpcUaParameters.RA)
        self.opc.get_positions.assert_called_once_with(OpcUaParameters.RA)

    def testMoveRelativeToCurrentPosition(self):
        self._mocked_pos = {}

        def side_effect_get(key):
            return self._mocked_pos[key]

        def side_effect_set(key, val):
            self._mocked_pos[key] = val

        self.opc.get_positions = Mock(side_effect=side_effect_get)
        self.opc.set_position = Mock(side_effect=side_effect_set)

        init_pos = 10
        self.angle.setPosition(init_pos)
        curr_pos = self.angle.getPosition()

        off = 19
        self.angle.setPosition(curr_pos + off)
        fin_pos = self.angle.getPosition()

        self.assertAlmostEqual(init_pos + off, fin_pos)

    def testOutOfRangeAngleRaises(self):
        with self.assertRaises(OSError):
            self.angle.setPosition(OpcUaParameters.max_angle + 1)
        with self.assertRaises(OSError):
            self.angle.setPosition(OpcUaParameters.min_angle - 1)
        self.opc.set_position.assert_not_called()
