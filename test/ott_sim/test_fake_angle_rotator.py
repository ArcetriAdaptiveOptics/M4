'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
from m4.ott_sim.fake_angle_rotator import FakeAngleRotator
#from m4.devices.base_angle_rotator import BaseAngleRotator


class TestFakeAngleRotator(unittest.TestCase):

    def testMoveRelativeToCurrentPosition(self):
        self.angle = FakeAngleRotator()
        init_pos = 42
        self.angle.setPosition(init_pos)

        off = 7
        curr_pos = self.angle.getPosition()
        fin_pos = self.angle.setPosition(curr_pos + off)

        self.assertAlmostEqual(init_pos + off, fin_pos)

#non riesco a testare il raise exception del BaseAngleRotator
#perchè non fa definire la classe senza le due definzioni:
#Can't instantiate abstract class AngleRotator with abstract methods getPosition, setPosition
