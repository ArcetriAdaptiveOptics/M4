import unittest
from unittest.mock import patch
import numpy as np
from arte.types.mask import CircularMask
from m4.type.measure import SurfaceMeasure
import datetime


# write test main for Measure class using mock on read_phasemap, _load_temperatures and _load_zernikes and fits.getheader
class testSurfaceMeasure(unittest.TestCase):

    def setUp(self):
        pass

    def test_load(self):

        with patch('m4.type.measure.read_phasemap') as mock_read_phasemap:
            mock_read_phasemap.return_value = np.ma.masked_array(
                np.ones((100, 100)), mask=CircularMask((100, 100), 50).mask)
            with patch('m4.type.measure.SurfaceMeasure._load_temperatures') as mock_load_temperatures:
                mock_load_temperatures.return_value = np.arange(
                    240).reshape(10, 24)
                with patch('m4.type.measure.SurfaceMeasure._load_zernikes') as mock_load_zernikes:
                    mock_load_zernikes.return_value = np.arange(
                        110).reshape(10, 11)
                    with patch('astropy.io.fits.getheader') as mock_getheader:
                        mock_getheader.return_value = {
                            'TN': '20240503_101010', 'N_MEAS': 4}
                        self.meas = SurfaceMeasure.load('20241231_235959.fits')

        self.assertEqual(self.meas.tn, '20241231_235959')
        self.assertEqual(self.meas.n_meas, 4)
        self.assertEqual(self.meas.temperatures.shape, (24,))
        self.assertEqual(self.meas.zernikes.shape, (11,))
        self.assertEqual(self.meas.shape, (100, 100))
        self.assertEqual(self.meas.timestamp, datetime.datetime.strptime(
            self.meas.tn, '%Y%m%d_%H%M%S'))
