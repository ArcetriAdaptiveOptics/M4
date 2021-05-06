'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import mock
from m4.devices.angle_rotator import OpcUaAngleRotator
from m4.configuration.ott_parameters import OpcUaParameters


class TestOpcUaAngleRotator(unittest.TestCase):

    def setUp(self):
        self.opc = mock.MagicMock()
        self.angle = OpcUaAngleRotator(self.opc)

    def testCallsOpcOnGetPosition(self):
        self.opc.get_position = mock.Mock(return_value=60)
        pos = self.angle.getPosition()
        self.opc.get_position.assert_called_once_with(OpcUaParameters.RA)
        self.assertAlmostEqual(pos, 60)

    def testCallsOpcOnSetPosition(self):
        self.opc.get_position = mock.Mock(return_value=3.14)
        pos = self.angle.setPosition(3.14)
        self.opc.set_target_position.assert_called_once_with(
            OpcUaParameters.RA, 3.14)
        self.opc.move_object.assert_called_once_with(OpcUaParameters.RA)
        self.opc.wait_for_stop.assert_called_once_with(OpcUaParameters.RA)
        self.assertAlmostEqual(pos, 3.14)

    def testMoveRelativeToCurrentPosition(self):

        self._mocked_pos = {}

        def side_effect_get(key):
            return self._mocked_pos[key]

        def side_effect_set(key, val):
            self._mocked_pos[key] = val
            print("%s" % self._mocked_pos)

        self.opc.get_position = mock.Mock(side_effect=side_effect_get)
        self.opc.set_target_position = mock.Mock(side_effect=side_effect_set)

        init_pos = 10
        self.angle.setPosition(init_pos)
        curr_pos = self.angle.getPosition()

        off = 19
        fin_pos = self.angle.setPosition(curr_pos + off)

        self.assertAlmostEqual(init_pos + off, fin_pos)
