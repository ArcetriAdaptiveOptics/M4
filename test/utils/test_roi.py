'''
Authors
  - C. Selmi: written in 2022
'''
import unittest
import os
from test.helper_test_library import testDataRootDir
from m4.ground import read_data

class TestRoi(unittest.TestCase):

    def _readImage(self, name):
        datadir = testDataRootDir()
        assert os.path.exists(datadir)
        directory = os.path.join(datadir, 'base',
                                              'M4Data', 'OPTData',
                                              'ROI')
        data_file_path = os.path.join(directory, name)
        image = read_data.readFits_maskedImage(data_file_path)
        return image

    def testROI(self):
        from m4.utils import roi
        r = roi.ROI()

        ima = self._readImage('seg_rmout.fits')
        dx, sx, c, rm = r.automatical_roi_selection(ima, True, False)

        ima = self._readImage('seg_rmin.fits')
        dx, sx, c, rm = r.automatical_roi_selection(ima, True, True)

        ima = self._readImage('cent_rmout.fits')
        segList = r.automatical_roi_selection(ima, False, False)

        ima = self._readImage('cent_rmin.fits')
        segList = r.automatical_roi_selection(ima, False, True)

        