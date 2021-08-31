'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
import mock
from m4.devices.reference_mirror_slider import OpcUaReferenceMirrorSlider
from m4.configuration.ott_parameters import OpcUaParameters

class TestOpcReferenceMirrorSlider(unittest.TestCase):

    def setUp(self):
        self.opc = mock.MagicMock()
        self.rslide = OpcUaReferenceMirrorSlider(self.opc)

    def testCallsOpcOnGetPosition(self):
        self.opc.get_position = mock.Mock(return_value=37)
        pos = self.rslide.getPosition()
        self.opc.get_position.assert_called_once_with(OpcUaParameters.CAR)
        self.assertAlmostEqual(pos, 37)

    def testCallsOpcOnSetPosition(self):
        self.opc.get_position = mock.Mock(return_value=100.14)
        pos = self.rslide.setPosition(100.14)
        self.opc.set_target_position.assert_called_once_with(
            OpcUaParameters.CAR, 100.14)
        self.opc.move_object.assert_called_once_with(OpcUaParameters.CAR)
        self.opc.wait_for_stop.assert_called_once_with(OpcUaParameters.CAR)
        self.assertAlmostEqual(pos, 100.14)

    def testMoveRelativeToCurrentPosition(self):
        self._mocked_pos = {}

        def side_effect_get(key):
            return self._mocked_pos[key]

        def side_effect_set(key, val):
            self._mocked_pos[key] = val
            print("%s" % self._mocked_pos)

        self.opc.get_position = mock.Mock(side_effect=side_effect_get)
        self.opc.set_target_position = mock.Mock(side_effect=side_effect_set)

        init_pos = 30
        self.rslide.setPosition(init_pos)
        curr_pos = self.rslide.getPosition()

        off = 19
        fin_pos = self.rslide.setPosition(curr_pos + off)

        self.assertAlmostEqual(init_pos + off, fin_pos)
