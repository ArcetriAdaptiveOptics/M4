'''
Authors
  - C. Selmi: written in 2021
'''
from unittest import TestCase
from unittest.mock import Mock, MagicMock
from m4.devices.reference_mirror_slider import OpcUaReferenceMirrorSlider
from m4.configuration.ott_parameters import OpcUaParameters

class TestOpcReferenceMirrorSlider(TestCase):

    def setUp(self):
        self.opc = MagicMock()
        self.rslide = OpcUaReferenceMirrorSlider(self.opc)

    def testCallsOpcOnGetPosition(self):
        self.opc.get_positions = Mock(return_value=37)
        pos = self.rslide.getPosition()
        self.opc.get_positions.assert_called_once_with(OpcUaParameters.CAR)
        self.assertAlmostEqual(pos, 37)

    def testCallsOpcOnSetPosition(self):
        self.opc.get_positions = Mock(return_value=100.14)
        # setPosition prints the reached position and returns nothing
        self.rslide.setPosition(100.14)
        self.opc.set_position.assert_called_once_with(OpcUaParameters.CAR, 100.14)
        self.opc.wait_for_stop.assert_called_once_with(OpcUaParameters.CAR)
        self.opc.get_positions.assert_called_once_with(OpcUaParameters.CAR)

    def testMoveRelativeToCurrentPosition(self):
        self._mocked_pos = {}

        def side_effect_get(key):
            return self._mocked_pos[key]

        def side_effect_set(key, val):
            self._mocked_pos[key] = val

        self.opc.get_positions = Mock(side_effect=side_effect_get)
        self.opc.set_position = Mock(side_effect=side_effect_set)

        init_pos = 30
        self.rslide.setPosition(init_pos)
        curr_pos = self.rslide.getPosition()

        off = 19
        self.rslide.setPosition(curr_pos + off)
        fin_pos = self.rslide.getPosition()

        self.assertAlmostEqual(init_pos + off, fin_pos)

    def testGetPositionInMeters(self):
        self.opc.get_positions = Mock(return_value=2500.0)
        self.assertAlmostEqual(self.rslide.getPositionInM(), 2.5)

    def testOutOfRangeSlideRaises(self):
        with self.assertRaises(ValueError):
            self.rslide.setPosition(OpcUaParameters.max_r_slide + 1)
        with self.assertRaises(ValueError):
            self.rslide.setPosition(OpcUaParameters.min_r_slide - 1)
        self.opc.set_position.assert_not_called()
