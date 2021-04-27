'''
Authors
  - C. Selmi: written in 2020
'''
import logging
import os
import time
import shutil
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits as pyfits
from m4.configuration.config import fold_name
from m4.devices.base_interferometer import BaseInterferometer

class I4dArcetri(BaseInterferometer):

    def __init__(self):
        """The constructor """
        from oaautils import i4d
        from m4.ground.read_data import InterferometerConverter
        self._ic = InterferometerConverter()
        self._interf = i4d.I4D()
        self._logger = logging.getLogger('4D')

    def acquire_phasemap(self, nframes=1, show=0):
        """
        Parameters
        ----------
            nframes: int
                number of frames

        Returns
        -------
            masked_ima: numpy masked array
                    interferometer image
        """
        if nframes == 1:
            masked_ima = self._getMeasurementOnTheFly(self._interf)
        else:
            cube_images = None
            for i in range(nframes):
                ima = self._getMeasurementOnTheFly(self._interf)
                if cube_images is None:
                    cube_images = ima
                else:
                    cube_images = np.ma.dstack((cube_images, ima))
            masked_ima = np.ma.mean(cube_images, axis=2)
        if show!=0:
            plt.clf()
            plt.imshow(masked_ima)
            plt.colorbar()
        return masked_ima

    def save_phasemap(self, dove, name, image):
        """
        Parameters
        ----------
        dove: string
            measurement file path
        name: string
            measuremnet fits file name
        image: numpy masked array
            data to save
        """
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, image.data)
        pyfits.append(fits_file_name, image.mask.astype(int))

    def _getMeasurementOnTheFly(self, interf):
        filename = '/tmp/prova4d'

        nMeasure=1
        interf.connect()
        interf.capture(1, name= 'DM_temp')
        interf.produce('DM_temp')
        interf.disconnect()
        time.sleep(1.0)
        fName = os.path.join(fold_name.PHASECAM_ROOT_FOLDER, 'DM_temp')
        #fName = '/home/m4/4d/Zcopy/DM_temp'

        for i in range(nMeasure):
            shutil.move(fName + '/hdf5/img_%04d.h5' %i,
                    filename + "_m%02d" %i + ".h5")

        shutil.rmtree(fName + '/hdf5')
        shutil.rmtree(fName + '/raw')

        return self._ic.from4D('/tmp/prova4d_m00.h5')
