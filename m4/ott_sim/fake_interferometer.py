'''
Authors
  - C. Selmi: written in 2020
'''
import logging
import os
import numpy as np
from astropy.io import fits as pyfits
from m4.ott_sim.ott_images import OttImages
from m4.devices.base_interferometer import BaseInterferometer

class FakeInterferometer(BaseInterferometer):

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('FakeInterferometer')
        self._ott = None

    def acquire_phasemap(self, n_frames=1, show=0):
        ottIma = OttImages(self._ott)
        opd, mask = ottIma.ott_smap(show=show)
        masked_ima = np.ma.masked_array(opd.T,
                                        mask=np.invert(mask.astype(bool)).T)
        return masked_ima

    def set_ott(self, ott):
        self._ott = ott

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

    def readImage4D(self, file_name):
        hduList = pyfits.open(file_name)
        masked_ima = np.ma.masked_array(hduList[0].data,
                                        hduList[1].data.astype(bool))
        return masked_ima
