'''
Authors
  - C. Selmi:  written in 2022
'''
import unittest
import os
import numpy as np
from m4.utils.influence_functions_maker import IFFunctionsMaker
from test.helper_test_library import testDataRootDir
from m4.configuration.start import create_ott
from m4.ott_sim.fake_interferometer import FakeInterferometer
from m4.utils.optical_calibration import OpticalCalibration
from m4.configuration.create_ott import OTT
import mock


class TestOpticalAlignment(unittest.TestCase):

    def setUp(self):
        self.cal = OpticalCalibration.loadCalibrationObjectFromFits(tt_cal)
        self.interf, self.ott = self._createInterferometer()

    def tetCreate(self):
        self.assertIsInstance(self.cal, OpticalCalibration)
        self.assertIsInstance(self.interf, FakeInterferometer)
        self.assertIsInstance(self.ott, OTT)


    @mock.patch('astropy.io.open', autospec=None)
    @mock.patch('m4.utils.influence_functions_maker.IFFunctionsMaker._storageFolder', autospec=True)
    @mock.patch('astropy.io.fits.writeto', autospec=None)
    @mock.patch('m4.type.modalAmplitude.ModalAmplitude._storageFolder', autospec=True)
    @mock.patch('m4.type.modalBase.ModalBase._storageFolder', autospec=True)
    @mock.patch('m4.ground.tracking_number_folder.os.makedirs', autospec=True)
    @mock.patch('m4.type.commandHistory.CmdHistory.saveInfo', autospec=True)

    def testReload(self, mock_folder):
        tt = '20220309_142454'
        mock_folder.return_value = os.path.join(testDataRootDir(), 'base',
                                                  'M4Data/OPTData/IFFunctions')
        IFFunctionsMaker.loadInfo(tt, 0)


    @mock.patch('m4.ground.read_data.readFits_data', autospec=True)
    @mock.patch('numpy.load', autospec=True)
    def _createInterferometer(self, mock_rd, mock_load):
        ott, interf, dm = create_ott(os.path.join(testDataRootDir(), 'base', 'Configurations', 'testConf.yaml'))
        return interf, ott
