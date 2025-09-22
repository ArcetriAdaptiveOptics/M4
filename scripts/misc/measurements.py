import time
import os
import shutil
import logging
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from opticalib.core import foldname as fname

class Measurements:
    """
    HOW TO USE IT::

    from m4.configuration import start
    ott, interf = start.create_ott(conf='.../youConf.yaml')
    from m4.mini_OTT.measurements import Measurements
    meas = Measurements(ott, interf)
    """

    def __init__(self, ott, interf):
        """The constructor"""
        self._ott = ott
        self._interf = interf


    def opticalMonitoring(self, n_images, delay, start_delay=0):
        """
        Acquisition of images for monitoring

        Parameters
        ----------
        n_images: int
            number of images to acquire
        delay: int [s]
            waiting time (in seconds) between two image acquisitions
        start_delay: int[s]
            waiting time before starting the measure (in seconds)

        Returns
        ------
        tt: string
            tracking number of measurements
        """
        store_in_folder = fold_name.OPD_SERIES_ROOT_FOLDER
        dove, tt = tracking_number_folder.createFolderToStoreMeasurements(
            store_in_folder
        )
        print(tt)
        shutil.copy(Interferometer.SETTINGS_CONF_FILE_M4OTT_PC, dove)
        shutil.move(
            os.path.join(dove, "AppSettings.ini"), os.path.join(dove, "4DSettings.ini")
        )
        ott_status.save_positions(dove, self._ott)  # saving the ott status

        print("waiting {:n} s...".format(start_delay))
        time.sleep(start_delay)
        print("start measuring")
        zer_list = []
        temp_list = []
        t0 = time.time()
        for i in range(n_images):
            ti = time.time()
            dt = ti - t0
            masked_ima = self._interf.acquire_phasemap(1)
            temp_vect = self._ott.temperature.getTemperature()
            name = Timestamp.now() + ".fits"
            fits_file_name = os.path.join(dove, name)
            pyfits.writeto(fits_file_name, masked_ima.data)
            pyfits.append(fits_file_name, masked_ima.mask.astype(np.uint8))

            coef, mat = zernike.zernikeFit(masked_ima, np.arange(10) + 1)
            vect = np.append(dt, coef)
            zer_list.append(vect)
            temp_list.append(temp_vect)

            fits_file_name = os.path.join(dove, "zernike.fits")
            pyfits.writeto(fits_file_name, np.array(zer_list), overwrite=True)
            fits_file_name = os.path.join(dove, "temperature.fits")
            pyfits.writeto(fits_file_name, np.array(temp_list), overwrite=True)

            time.sleep(delay)
        return tt

