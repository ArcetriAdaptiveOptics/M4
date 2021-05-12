'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
from m4.ott_sim.fake_angle_rotator import FakeAngleRotator


class TestFakeAngleRotator(unittest.TestCase):

    def testMoveRelativeToCurrentPosition(self):
        self.angle = FakeAngleRotator()
        init_pos = 42
        self.angle.setPosition(init_pos)

        off = 7
        curr_pos = self.angle.getPosition()
        fin_pos = self.angle.setPosition(curr_pos + off)

        self.assertAlmostEqual(init_pos + off, fin_pos)
