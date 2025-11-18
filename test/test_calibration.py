'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import numpy as np
import os
from unittest.mock import patch
from test.helper_test_library import testDataRootDir
# Note: m4.utils.optical_calibration.OpticalCalibration has been removed.
# The functionality has been moved to opticalib.alignment.Alignment.
# These tests are kept for reference but need to be rewritten to use the new API.
try:
    from m4.utils import optical_calibration
except ImportError:
    optical_calibration = None

PARDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'img_0000.fits')

class TestCalc(unittest.TestCase):

    def setUp(self):
        if optical_calibration is None:
            self.skipTest('m4.utils.optical_calibration has been removed. '
                         'Functionality moved to opticalib.alignment.Alignment. '
                         'Tests need to be updated to use the new API.')
        self.ott, self.interf = self._createOttAndInterf()
        self.cal = optical_calibration.OpticalCalibration(self.ott, self.interf)

    def tearDown(self):
        del(self.cal, self.ott)

    @unittest.skip('Da rifare')
    @patch('m4.ground.read_data.readFits_data', autospec=True)
    @patch('numpy.load', autospec=True)
    def _createOttAndInterf(self, mock_rd, mock_load):
        from m4.configuration.ott import create_ott
        from m4.simulator.fake_parabola_slider import FakeParabolaSlider
        from m4.simulator.fake_interferometer import FakeInterferometer
        ott, interf, dm = create_ott()
        self.assertIsInstance(ott.parabolaSlider, FakeParabolaSlider)
        self.assertIsInstance(interf, FakeInterferometer)
        interf.save_phasemap = self._skipSave
        interf.acquire_phasemap = self._skipAcq
        
        return ott, interf

    def _skipAcq(self, aa, bb):
        image = np.zeros((500, 500))
        mask = np.ones((500, 500), dtype=bool)
        masked_ima = np.ma.masked_array(image, mask=mask)
        return masked_ima

    def _skipSave(self, a, b, c):
        pass

    def _skipParRemapped(self, a, b):
        image = np.zeros((500, 500))
        mask = np.ones((500, 500), dtype=bool)
        masked_ima = np.ma.masked_array(image, mask=mask)
        return masked_ima

    @unittest.skipIf(optical_calibration is None, 
                     'm4.utils.optical_calibration has been removed. '
                     'Functionality moved to opticalib.alignment.Alignment.')
    @patch('astropy.io.fits.writeto', autospec=None)
    @patch('m4.ground.tracking_number_folder._error', autospec=None)
    @patch('m4.ground.timestamp.Timestamp.now', autospec=True)
    @patch('m4.utils.optical_calibration.OpticalCalibration._storageFolder', autospec=True)
    @patch('m4.ground.tracking_number_folder.os.makedirs', autospec=True)
    def testCalibration(self, mock_makedirs, mockFilepath1, mock_tt, mock_oserror, mock_savecalib):
        who = optical_calibration.WHO_PAR_AND_RM
        command_amp_vector = np.ones(5)
        n_push_pull = 1

        mock_tt.return_value = '20220217_151631'
        mockFilepath1.return_value = os.path.join(testDataRootDir(), 'base',
                                                  'M4Data/OPTData/Calibration')
        self.cal.measureAndAnalysisCalibrationMatrix(who, command_amp_vector,
                                                     n_push_pull, n_frames=None, delay=None)

    @unittest.skipIf(optical_calibration is None, 
                     'm4.utils.optical_calibration has been removed. '
                     'Functionality moved to opticalib.alignment.Alignment.')
    @patch('m4.utils.optical_calibration.OpticalCalibration._storageFolder', autospec=True)
    def testReload(self, mockFilepath1):
        mockFilepath1.return_value = os.path.join(testDataRootDir(), 'base',
                                                  'M4Data/OPTData/Calibration')
        tt = '20220217_151631'
        cal = optical_calibration.OpticalCalibration.loadCalibrationObjectFromFits(tt)
        cal.getCube()
        cal.getMask()
        cal.getWho()
