"""
Authors
  - C. Selmi:  written in 2022
"""

import unittest
import os
import numpy as np
from m4.utils.influence_functions_maker import IFFunctionsMaker
from m4.configuration.start import create_ott
from m4.ott_sim.fake_interferometer import FakeInterferometer
from m4.utils.optical_calibration import OpticalCalibration
from m4.configuration.create_ott import OTT as ottobj
from m4.ott_sim.fake_deformable_mirror import FakeM4DM
from m4.utils.optical_alignment import OpticalAlignment
import mock
from unittest.mock import MagicMock


class TestOpticalAlignment(unittest.TestCase):

    class FakeOpen:
        def __init__(self, filename, mode):
            self.filename = filename
            self.mode = mode
            self.write_calls = 0

        def read(self):
            return b"{}"

        def write(self, write_string):
            self.write_calls += 1

        def close(self):
            if self.mode == "w":
                assert self.write_calls == 1

        def __enter__(self):
            return self

        def __exit__(self, *args, **kwargs):
            pass

    def testDataRootDir(self):
        return os.path.join(os.path.dirname(__file__), "../data")

    def tearDown(self):
        self._cal
        self._ott
        self._interf
        self._dm

    @mock.patch("m4.utils.optical_calibration.OpticalCalibration", autospec=True)
    @mock.patch("astropy.io.fits.getheader", autospec=True)
    @mock.patch("astropy.io.fits.open", autospec=True)
    @mock.patch("m4.configuration.start.FakeM4DM", autospec=True)
    def setUp(self, mock_interf, mock_open, mock_header, mock_cal):
        self._fake_header = {"WHO": "M4", "NPUSHPUL": 1}
        self._fake_data = [np.ones((10, 10)) for i in range(4)]
        tt_cal = "20220309_142454.fits"
        mock_header.return_value = self._fake_header

        # Create a mock HDUList object
        mock_hdulist = MagicMock()
        mock_hdulist.__enter__.return_value = mock_hdulist
        mock_hdulist.__exit__.return_value = False
        mock_hdulist[0].data = np.ones((10, 10))
        mock_open.return_value = mock_hdulist

        mock_cal.return_value = OpticalCalibration(mock_interf, None)
        self._cal = OpticalCalibration.loadCalibrationObjectFromFits(tt_cal)
        self.assertIsInstance(self._cal, OpticalCalibration)

        self._ott, self._interf, self._dm = create_ott()
        self.assertIsInstance(self._ott, ottobj)
        self.assertIsInstance(self._interf, FakeInterferometer)
        self.assertIsInstance(self._dm, FakeM4DM)

    @mock.patch("m4.utils.optical_calibration.OpticalCalibration", autospec=True)
    @mock.patch("astropy.io.fits.getheader", autospec=True)
    @mock.patch("astropy.io.fits.open", autospec=True)
    @mock.patch.object(
        OpticalAlignment,
        "selectModesInIntMatAndRecConstruction",
        return_value=(np.identity(10), np.identity(10), np.identity(10)),
    )
    @mock.patch("os.makedirs", autospec=True)
    @mock.patch("m4.utils.optical_alignment.open", FakeOpen, create=True)
    @mock.patch("m4.configuration.start.FakeM4DM", autospec=True)
    @mock.patch("m4.configuration.start.create_ott", autospec=True)
    @mock.patch.object(
        FakeInterferometer,
        "acquire_phasemap",
        return_value=np.ma.masked_array(np.ones((10, 10))),
    )
    @mock.patch.object(
        FakeInterferometer,
        "save_phasemap",
        return_value=np.ma.masked_array(np.ones((10, 10))),
    )
    def testOptAlign(
        self,
        mock_sphm,
        mock_aphm,
        mock_cott,
        mock_interf,
        mock_makedirs,
        mock_selectModes,
        mock_fits_open,
        mock_header,
        mock_cal,
    ):
        self._ott, self._interf, self._dm = create_ott()
        tt_cal = "20220309_142454.fits"
        self._align = OpticalAlignment(tt_cal, self._ott, FakeInterferometer(self._ott, self._dm))
        self._align.opt_aligner(
            n_images=10,
            delay=1,
            zernike_to_be_corrected=[1, 2, 3],
            dof_command_id=None,
            subapOffsets=False,
        )
