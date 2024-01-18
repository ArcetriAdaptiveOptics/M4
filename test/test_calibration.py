'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import numpy as np
import os
import mock
from test.test_helper import testDataRootDir
from m4.utils import optical_calibration

PARDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'img_0000.fits')

class TestCalc(unittest.TestCase):

    def setUp(self):
        self.ott, self.interf = self._createOttAndInterf()
        self.cal = optical_calibration.OpticalCalibration(self.ott, self.interf)

    def tearDown(self):
        del(self.cal, self.ott)

    @unittest.skip('Da rifare')
    @mock.patch('m4.ground.read_data.readFits_data', autospec=True)
    @mock.patch('numpy.load', autospec=True)
    def _createOttAndInterf(self, mock_rd, mock_load, mockFilepath):
        from m4.configuration.start import create_ott
        from m4.ott_sim.fake_parabola_slider import FakeParabolaSlider
        from m4.ott_sim.fake_interferometer import FakeInterferometer
        ott, interf, dm = create_ott(os.path.join(testDataRootDir(), 'base',
                                              'Configurations', 'testConf.yaml'))
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

    @mock.patch('astropy.io.fits.writeto', autospec=None)
    @mock.patch('m4.ground.tracking_number_folder._error', autospec=None)
    @mock.patch('m4.ground.timestamp.Timestamp.now', autospec=True)
    @mock.patch('m4.utils.optical_calibration.OpticalCalibration._storageFolder', autospec=True)
    @mock.patch('m4.ground.tracking_number_folder.os.makedirs', autospec=True)
    def testCalibration(self, mock_makedirs, mockFilepath1, mock_tt, mock_oserror, mock_savecalib):
        who = optical_calibration.WHO_PAR_AND_RM
        command_amp_vector = np.ones(5)
        n_push_pull = 1

        mock_tt.return_value = '20220217_151631'
        mockFilepath1.return_value = os.path.join(testDataRootDir(), 'base',
                                                  'M4Data/OPTData/Calibration')
        self.cal.measureAndAnalysisCalibrationMatrix(who, command_amp_vector,
                                                     n_push_pull, n_frames=None, delay=None)

    @mock.patch('m4.utils.optical_calibration.OpticalCalibration._storageFolder', autospec=True)
    def testReload(self, mockFilepath1):
        mockFilepath1.return_value = os.path.join(testDataRootDir(), 'base',
                                                  'M4Data/OPTData/Calibration')
        tt = '20220217_151631'
        cal = optical_calibration.OpticalCalibration.loadCalibrationObjectFromFits(tt)
        cal.getCube()
        cal.getMask()
        cal.getWho()
