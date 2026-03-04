import time
import os
import shutil
import logging
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from opticalib.ground.osutils import load_fits, save_fits, newtn
from opticalib.ground import modal_decomposer as zern
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

    def __init__(self, interf, ott=None, surfFitRoutine = None):
        """The constructor"""
        self._ott = ott
        self._interf = interf
        self._sfit = surfFitRoutine  #così??
        self.basepath = fold_name.BASE_DATA_PATH


    def opticalMonitoring(self, n_images, delay=0, start_delay=0, fullFrame=False, tracknum = None):
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
        if tracknum is None:
            tt = newtn()
        else:
            tt = tracknum
        savefolder = os.path.join(store_in_folder, tt)
        print(tt)
        if os.path.exists(savefolder) == False:
            os.mkdir(savefolder)
        self._interf.copy4DSettings(savefolder)
        #shutil.copy(Interferometer.SETTINGS_CONF_FILE_M4OTT_PC, dove)
        #shutil.move(os.path.join(dove, "AppSettings.ini"), os.path.join(dove, "4DSettings.ini")
        #if self._ott is not None:
        #    ott_status.save_positions(savefolder, self._ott)  # saving the ott status

        print("waiting {:n} s...".format(start_delay))
        time.sleep(start_delay)
        print("start measuring")
        zer_list = []
        temp_list = []
        t0 = time.time()
        for i in range(n_images):
            ti = time.time()
            dt = ti - t0
            masked_ima = self._interf.acquire_map(1)
            if fullFrame is not False:
                masked_ima = self._interf.intoFullFrame(1)
            if self._ott is not None:
                temp_vect = self._ott.temperature.getTemperature()
                fits_file_name = os.path.join(savefolder, "temperature.fits")
                opt.save_fits(fits_file_name, np.array(temp_list), overwrite=True)
                temp_list.append(temp_vect)

            name = newtn() + ".fits"
            fits_file_name = os.path.join(savefolder, name)
            opt.save_fits(fits_file_name, masked_ima)
            if self._sfit is not None:
                sfits_file_name = os.path.join(savefolder, "zernike.fits")
                coeff = self._sfit(masked_ima)
                coeff_list.append(dt,coeff)
                opt.save_fits(sfits_file_name, np.array(coeff_list), overwrite=True)

            #coef, mat = zern.zernikeFit(masked_ima, np.arange(10) + 1) ## qui bisognerebbe capire come impostare il fit di zernike
            #vect = np.append(dt, coef)
            #zer_list.append(vect)

            time.sleep(delay)
        return tt

