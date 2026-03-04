import time
import os
import shutil
import logging
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
#from opticalib import measurements as measlib
from scripts.misc import measurements as measlib
from opticalib.ground.osutils import load_fits, save_fits, newtn
from opticalib.ground import modal_decomposer as zern
import opticalib as opt
from m4.configuration import ott_status

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
        self.meas = measlib.Measurements(interf, ott)

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
        #self._interf.copy4DSettings(savefolder)   #in meas
        #shutil.copy(Interferometer.SETTINGS_CONF_FILE_M4OTT_PC, dove)
        #shutil.move(os.path.join(dove, "AppSettings.ini"), os.path.join(dove, "4DSettings.ini")
        if self._ott is not None:
            ott_status.save_positions(savefolder, self._ott)  # saving the ott status
        self.meas.opticalMonitoring(n_images, delay, start_delay, fullFrame, tt)
        return tt

