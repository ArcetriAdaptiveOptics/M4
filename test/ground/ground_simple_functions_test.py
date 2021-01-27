'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import numpy as np
from m4.ground import geo
from m4.ground import find_directory
from m4.ground import smooth_function
from m4.ground.timestamp import Timestamp


class Test(unittest.TestCase):


    @unittest.skip('Mettere dei file')
    def testFindDirectory(self):
        tt = '2021...'
        find_directory.findTtPath(tt)

    def testGeometry(self):
        img = np.random.rand(500, 500)
        masked_ima = geo.draw_mask(img, 250, 250, 50)
        geo.qpupil(masked_ima)
        geo.rotate(img, 30)

    def testSmoothFunction(self):
        data = np.arange(100)
        smooth_function.smooth(data, 4)

    def testTimestamp(self):
        t = Timestamp()
        t.asNowString()
        t.asTodayString()
        t.now()
        t.nowUSec()
        t.today()