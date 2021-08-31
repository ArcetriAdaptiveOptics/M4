'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
import numpy as np
from numpy import testing
from m4.ott_sim.fake_parabola import FakeParabola

class TestFakeReferenceMirror(unittest.TestCase):

    def testMoveRelativeToCurrentPosition(self):
        self.parabola = FakeParabola()
        init_pos = np.array([0, 0, 3e-3, 1e-6, 2e-7, 0])
        self.parabola.setPosition(init_pos)

        off = np.array([0, 0, 1e-3, 3e-6, -7e-7, 0])
        curr_pos = self.parabola.getPosition()
        fin_pos = self.parabola.setPosition(curr_pos + off)

        testing.assert_almost_equal(init_pos + off, fin_pos)
