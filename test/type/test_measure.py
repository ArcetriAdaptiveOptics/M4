import unittest
from unittest.mock import patch
import numpy as np
from arte.types.mask import CircularMask
from m4.type.measure import SurfaceMeasure


# write test main for Measure class using mock on read_phasemap, _load_temperatures and _load_zernikes and fits.getheader
class testSurfaceMeasure(unittest.TestCase):

    def setUp(self):
        pass

    def test_load(self):

        with patch('m4.type.measure.read_phasemap') as mock_read_phasemap:
            mock_read_phasemap.return_value = np.ma.masked_array(
                np.ones((100, 100)), mask=CircularMask((100, 100), 50).mask)
            with patch('m4.type.measure.SurfaceMeasure._load_temperatures') as mock_load_temperatures:
                mock_load_temperatures.return_value = np.arange(240).reshape(24, 10)
                with patch('m4.type.measure.SurfaceMeasure._load_zernikes') as mock_load_zernikes:
                    mock_load_zernikes.return_value = np.arange(110).reshape(11, 10)
                    with patch('astropy.io.fits.getheader') as mock_getheader:
                        mock_getheader.return_value = {
                            'TN': '20240503_101010', 'N_MEAS': 4}
                        meas = SurfaceMeasure.load('20230715_151001.fits')
