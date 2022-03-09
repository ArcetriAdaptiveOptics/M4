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

class TestCalc(unittest.TestCase):

    def setUp(self):
        self.ott, self.interf = self._createOttAndInterf()
        self.cal = optical_calibration.OpticalCalibration(self.ott, self.interf)

    def tearDown(self):
        del(self.cal, self.ott)

    def _createOttAndInterf(self):
        from m4.configuration.start import create_ott
        from m4.ott_sim.fake_parabola_slider import FakeParabolaSlider
        from m4.ott_sim.fake_interferometer import FakeInterferometer
        ott, interf = create_ott(os.path.join(testDataRootDir(), 'base',
                                              'Configurations', 'testConf.yaml'))
        self.assertIsInstance(ott.parabolaSlider, FakeParabolaSlider)
        self.assertIsInstance(interf, FakeInterferometer)
        interf.save_phasemap = self._skipSave
        interf.acquire_phasemap = self._skipAcq
        return ott, interf

    def _skipAcq(self, aa):
        pass

    def _skipSave(self, a, b, c):
        pass


    @mock.patch('astropy.io.fits.writeto', autospec=None)
    @mock.patch('m4.ground.tracking_number_folder._error', autospec=None)
    @mock.patch('m4.ground.timestamp.Timestamp.now', autospec=True)
    @mock.patch('m4.utils.optical_calibration.OpticalCalibration._storageFolder', autospec=True)
    @mock.patch('m4.ground.tracking_number_folder.os.makedirs', autospec=True)
    def testCalibration(self, mock_makedirs, mockFilepath1, mock_tt, mock_oserror, mock_savecalib):
        who = optical_calibration.WHO_PAR_AND_RM
        command_amp_vector = np.ones(5)
        n_push_pull = 1
        n_frames = 1

        mock_tt.return_value = '20220217_151631'
        mockFilepath1.return_value = os.path.join(testDataRootDir(), 'base',
                                                  'M4Data/OPTData/Calibration')
        self.cal.measureAndAnalysisCalibrationMatrix(who, command_amp_vector,
                                                     n_push_pull, n_frames)

    @mock.patch('m4.utils.optical_calibration.OpticalCalibration._storageFolder', autospec=True)
    def testReload(self, mockFilepath1):
        mockFilepath1.return_value = os.path.join(testDataRootDir(), 'base',
                                                  'M4Data/OPTData/Calibration')
        tt = '20220217_151631'
        cal = optical_calibration.OpticalCalibration.loadCalibrationObjectFromFits(tt)
        cal.getCube()
        cal.getMask()
        cal.getWho()
