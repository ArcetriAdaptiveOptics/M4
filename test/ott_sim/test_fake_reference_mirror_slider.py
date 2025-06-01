'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
from m4.simulator.fake_reference_mirror_slider import FakeReferenceMirrorSlider


class TestFakeReferenceMirrorSlider(unittest.TestCase):

    def testMoveRelativeToCurrentPosition(self):
        self.rslide = FakeReferenceMirrorSlider()
        init_pos = 77
        self.rslide.setPosition(init_pos)

        off = 3
        curr_pos = self.rslide.getPosition()
        fin_pos = self.rslide.setPosition(curr_pos + off)

        self.assertAlmostEqual(init_pos + off, fin_pos)
