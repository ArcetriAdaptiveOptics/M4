'''
Authors
  - C. Selmi:  written in 2022
'''
import unittest
import os
from m4.analyzers import requirement_analyzer as req_check
from test.helper_test_library import testDataRootDir
from m4.ground import read_data

class TestRequirementAnalyzer(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @unittest.skip('Da rifare')
    def testReadImageAndAnalysis(self):
        image_location = os.path.join(testDataRootDir(),
                                      'ima_for_req.fits')
        image = read_data.read_phasemap(image_location)

        pscale = None
        step = 10000
        n_patches = None
        slope = req_check.test242(image, pscale)
        diff_piston = req_check.diffPiston(image)
        roc = req_check.test283(image, pscale, step)
        rem_31 = req_check.test243(image, 0.015, pscale, step, n_patches)
        rms_500 = req_check.test243(image, 0.1, pscale, step, n_patches)
        