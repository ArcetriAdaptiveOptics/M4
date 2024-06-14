"""
Authors
  - C. Selmi:  written in 2022
"""

import unittest
import os
import numpy as np
from m4.utils.influence_functions_maker import IFFunctionsMaker
from test.helper_test_library import testDataRootDir
from m4.configuration.start import create_ott
from m4.ott_sim.fake_interferometer import FakeInterferometer
from m4.utils.optical_calibration import OpticalCalibration
from m4.configuration.create_ott import OTT as ottobj
import mock


class TestOpticalAlignment(unittest.TestCase):


    def optical_calibration_init_mock(self, tt, dm, interf, ott):
        """The mock  constructor """
        self._logger = logging.getLogger('OPT_CALIB:')
        self._interf = interf
        self._ott = ott
        #start
        self._nPushPull = None
        self._commandAmpVector = None
        self._who = None
        #from calibration
        self._dofIndex = None
        self._commandMatrix = None
        self._commandList = None
        self.tt = None
        #from analysis
        self._cube = None
        self._mask = None
        self._intMat = None

        self._fullCommandMatrix=None
        self._fullCube=None

    def setUp(self):
        self._fake_header = {"WHO": "M4", "NPUSHPUL": 1}
        self._fake_data = [np.ones((10, 10)) for i in range(4)]

    @mock.patch.object(m4.utils.optical_calibration.OpticalCalibration, "__init__", optical_calibration_init_mock)
    @mock.patch("m4.utils.optical_calibration.OpticalCalibration",autospec=True)
#    @mock.patch("m4.configuration.start", autospec=True)
    @mock.patch("astropy.io.fits.getheader", autospec=True)
    @mock.patch("astropy.io.fits.open", autospec=True)
    def testCreate(self, mock_cal, mock_header, mock_open):
        tt_cal = "20220309_142454.fits"
        mock_header.return_value = self._fake_header
        mock_open.return_value = self._fake_data
        self.cal = OpticalCalibration.loadCalibrationObjectFromFits(tt_cal)
        self.interf, self.ott, self.dm = create_ott(
            os.path.join(testDataRootDir(), "base", "Configurations", "testConf.yaml")
        )

    #    self.assertIsInstance(self.interf, FakeInterferometer)
    #    self.assertIsInstance(self.ott, ottobj)
    #    self.dom(self.dm, dm)
    #    pass


#    @mock.patch(
#        "m4.utils.influence_functions_maker.IFFunctionsMaker._storageFolder",
#        autospec=True,
#    )
#    @mock.patch("astropy.io.fits.writeto", autospec=None)
#    @mock.patch("m4.type.modalAmplitude.ModalAmplitude._storageFolder", autospec=True)
#    @mock.patch("m4.type.modalBase.ModalBase._storageFolder", autospec=True)
#    @mock.patch("m4.ground.tracking_number_folder.os.makedirs", autospec=True)
#    @mock.patch("m4.type.commandHistory.CmdHistory.saveInfo", autospec=True)
#    def testReload(self, mock_folder):
#        tt = "20220309_142454"
#        mock_folder.return_value = os.path.join(
#            testDataRootDir(), "base", "M4Data/OPTData/IFFunctions"
#        )
#        IFFunctionsMaker.loadInfo(tt, 0)

#write main class to run setUp and testCreate
if __name__ == "__main__":
    unittest.main()

