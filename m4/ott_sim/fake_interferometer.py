'''
Authors
  - C. Selmi: written in 2020
'''
import logging
import os
import time
import numpy as np
from astropy.io import fits as pyfits
from m4.ott_sim.ott_images import OttImages
from m4.devices.base_interferometer import BaseInterferometer

class FakeInterferometer(BaseInterferometer):
    '''

    HOW TO USE IT::

        from m4.ott_sim.fake_interferometer import FakeInterferometer
        interf = FakeInterferometer()
        image = interf.acquire_phasemap()
    '''

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('FakeInterferometer')
        self._ott = None
        self._dm = None

    def acquire_phasemap(self, n_frames=1, delay=0, indet=True):
        '''
        Parameters
        ----------
            nframes: int
                number of frames
            delay: int [s]
                delay between images

        Other Parameters
        ----------------
            indet: boolean
                True to consider lamba indeterminacy

        Returns
        -------
            masked_ima: numpy masked array
                    interferometer image
        '''
        ima_list = []
        for i in range(n_frames):
            ottIma = OttImages(self._ott)
            ottIma.m4ima = self._dm.m4ima
            opd, mask = ottIma.ott_smap(show=0)
            masked_ima = np.ma.masked_array(opd.T,
                                            mask=np.invert(mask.astype(bool)).T)
            '''
            aggiungo indeterminazione di lambda
            Luca
            '''
            if indet==True:
                lam = 632.8e-9
                kk = np.floor(np.random.random(1)*5-2) 
                masked_ima = masked_ima + np.ones(masked_ima.shape)*lam*kk
            ima_list.append(masked_ima)
            time.sleep(delay)
        images = np.dstack(ima_list)
        ima = np.mean(images, 2)
        masked_ima = np.ma.masked_array(ima, mask=np.invert(mask.astype(bool)).T)

        return masked_ima

    def set_ott(self, ott):
        ''' Function for setting optical tower data
        '''
        self._ott = ott

    def set_dm(self, deformable_mirror):
        ''' Function for setting deformable mirror data
        '''
        self._dm = deformable_mirror

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
