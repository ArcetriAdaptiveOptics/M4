import time
import os
import shutil
import logging
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from opticalib import measurements as meas
from opticalib.ground.osutils import load_fits, save_fits, newtn
from opticalib.ground import zernike as zern
import opticalib as opt

fold_name = opt.core.root.folders

class Measurements:
    """
    HOW TO USE IT::

    from m4.configuration import start
    ott, interf = start.create_ott(conf='.../youConf.yaml')
    from m4.mini_OTT.measurements import Measurements
    meas = Measurements(ott, interf)
    """

    def __init__(self, interf, ott=None):
        """The constructor"""
        self._ott = ott
        self._interf = interf
        self.basepath = fold_name.BASE_DATA_PATH
        self.meas = meas(interf, ott)

    def opticalMonitoring(self, n_images, delay=0, start_delay=0, fullFrame=False):
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
        tt = newtn()
        savefolder = os.path.join(store_in_folder, tt)
        print(tt)
        if os.path.exists(savefolder) == False:
            os.mkdir(savefolder)
        #self._interf.copy4DSettings(savefolder)   #in meas
        #shutil.copy(Interferometer.SETTINGS_CONF_FILE_M4OTT_PC, dove)
        #shutil.move(os.path.join(dove, "AppSettings.ini"), os.path.join(dove, "4DSettings.ini")
        if self._ott is not None:
            ott_status.save_positions(savefolder, self._ott)  # saving the ott status
        self.meas(opticalMonitoring(n_images, delay, start_delay, fullFrame, tt)
        '''
        print("waiting {:n} s...".format(start_delay))
        time.sleep(start_delay)
        print("start measuring")
        zer_list = []
        temp_list = []
        t0 = time.time()
        for i in range(n_images):
            ti = time.time()
            dt = ti - t0
#            masked_ima = self._interf.acquire_phasemap(1)
            masked_ima = self._interf.acquire_map(1)
            if fullFrame is not False:
                masked_ima = self._interf.intoFullFrame(1)
            if self._ott is not None:
                temp_vect = self._ott.temperature.getTemperature()
                fits_file_name = os.path.join(savefolder, "temperature.fits")
                opt.save_fits(fits_file_name, np.array(temp_list), overwrite=True)

            name = newtn() + ".fits"
            fits_file_name = os.path.join(savefolder, name)
            opt.save_fits(fits_file_name, masked_ima)
 #           pyfits.writeto(fits_file_name, masked_ima.data)
 #           pyfits.append(fits_file_name, masked_ima.mask.astype(np.uint8))

            coef, mat = zernike.zernikeFit(masked_ima, np.arange(10) + 1)
            vect = np.append(dt, coef)
            zer_list.append(vect)
            temp_list.append(temp_vect)

            fits_file_name = os.path.join(savefolder, "zernike.fits")
            opt.save_fits(fits_file_name, np.array(zer_list), overwrite=True)
#            pyfits.writeto(fits_file_name, np.array(zer_list), overwrite=True)
#            pyfits.writeto(fits_file_name, np.array(temp_list), overwrite=True)

            time.sleep(delay)
            '''
        return tt

