'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
import numpy as np
from numpy import testing
from m4.ott_sim.fake_reference_mirror import FakeReferenceMirror

class TestFakeReferenceMirror(unittest.TestCase):

    def testMoveRelativeToCurrentPosition(self):
        self.rm = FakeReferenceMirror()
        init_pos = np.array([0, 0, 1e-13, 3e-6, 7e-7, 0])
        self.rm.setPosition(init_pos)

        off = np.array([0, 0, 0, 3e-6, -7e-7, 0])
        curr_pos = self.rm.getPosition()
        fin_pos = self.rm.setPosition(curr_pos + off)

        testing.assert_almost_equal(init_pos + off, fin_pos)
