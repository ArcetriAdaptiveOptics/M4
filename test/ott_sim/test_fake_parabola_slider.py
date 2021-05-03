#!/usr/bin/env python
import unittest
from m4.ott_sim.fake_parabola_slider import FakeParabolaSlider


class TestFakeParabolaSlider(unittest.TestCase):

    def testMoveRelativeToCurrentPosition(self):
        self.slider = FakeParabolaSlider()
        init_pos = 42
        self.slider.setPosition(init_pos)

        off = 19.3
        curr_pos = self.slider.getPosition()
        fin_pos = self.slider.setPosition(curr_pos + off)

        self.assertAlmostEqual(init_pos + off, fin_pos)
