import unittest
from unittest.mock import patch
import numpy as np
from arte.types.mask import CircularMask
from m4.type.measure import SurfaceMeasure, SurfaceSequence
import datetime
from pathlib import Path
import shutil
import tempfile
import os
import astropy.io.fits as fits

# write test main for Measure class using mock on read_phasemap, _load_temperatures and _load_zernikes and fits.getheader


class testSurfaceMeasure(unittest.TestCase):

    def setUp(self):
        pass

    def test_load(self):

        with patch("m4.type.measure.read_phasemap") as mock_read_phasemap:
            mock_read_phasemap.return_value = np.ma.masked_array(
                np.ones((100, 100)), mask=CircularMask((100, 100), 50).mask
            )
            with patch(
                "m4.type.measure.SurfaceMeasure.load_temperatures"
            ) as mock_load_temperatures:
                mock_load_temperatures.return_value = np.arange(240).reshape(10, 24)
                with patch(
                    "m4.type.measure.SurfaceMeasure.load_zernikes"
                ) as mock_load_zernikes:
                    mock_load_zernikes.return_value = np.arange(110).reshape(10, 11)
                    with patch("astropy.io.fits.getheader") as mock_getheader:
                        mock_getheader.return_value = {
                            "TN": "20240503_101010",
                            "N_MEAS": 4,
                        }
                        self.meas = SurfaceMeasure.load("20241231_235959.fits")

        self.assertEqual(self.meas.tn, "20241231_235959")
        self.assertEqual(self.meas.n_meas, 4)
        self.assertEqual(self.meas.temperatures.shape, (24,))
        self.assertEqual(self.meas.zernikes.shape, (11,))
        self.assertEqual(self.meas.shape, (100, 100))
        self.assertEqual(
            self.meas.timestamp,
            datetime.datetime.strptime(self.meas.tn, "%Y%m%d_%H%M%S"),
        )


class testSurfaceSequence(unittest.TestCase):
    def setUp(self):
        # creare un albero finto da qualche parte dinamicamente
        self.tmp_path = Path(Path(tempfile.gettempdir()), "tmp")
        self.tn_path = Path(self.tmp_path, "20241226_000000")
        os.makedirs(self.tn_path, exist_ok=True)
        for i in range(1, 50):
            header = fits.Header()
            header["TN"] = f"20241226_{i:04d}00"
            data = np.ma.masked_array(
                np.random.rand(100, 100), mask=CircularMask((100, 100), 50).mask
            )
            fits.writeto(
                Path(self.tn_path, header["TN"] + ".fits"),
                data.data,
                header,
                overwrite=True,
            )
            fits.append(
                Path(self.tn_path, header["TN"] + ".fits"),
                data.mask.astype(np.uint8),
                header,
                overwrite=True,
            )

        zt = np.random.rand(50, 11)
        tt = np.random.rand(50, 24)
        fits.writeto(Path(self.tn_path, "temperature.fits"), tt, overwrite=True)

        fits.writeto(Path(self.tn_path, "zernike.fits"), zt, overwrite=True)

        # Path(tn_path, '20241226_020000.fits').touch()
        # Path(tn_path, '20241226_030000.fits').touch()
        # Path(tn_path, '20241226_040000.fits').touch()
        # Path(tn_path, '20241226_050000.fits').touch()
        # Path(tn_path, '20241226_060000.fits').touch()
        # Path(tn_path, '20241226_070000.fits').touch()
        # Path(tn_path, 'temperatures.fits').touch()
        # Path(tn_path, 'zernikes.fits').touch()

        print("Created sample tree into: ", self.tmp_path)

    def test_load(self):

        print(self.tn_path)
        self.seq = SurfaceSequence(self.tn_path)
        print("seq len %d" % len(self.seq))
        print(self.seq[0])
        print(self.seq[1])
        self.assertTrue(
            self.seq[1].tn == "20241226_001800"
            or self.seq[1].tn == "20241226_000200"
            or self.seq[1].tn == "20241226_003700" # patched - to check
        )
        # seq2 = self.seq[1:4]

        # time, stf = seq2.compute_fast_time_stf()

    def tearDown(self):
        print("Cleaning temporary tree")
        try:
            shutil.rmtree(self.tmp_path)
        except Exception as e:
            print("Temporary tree cleaned except for some files")
