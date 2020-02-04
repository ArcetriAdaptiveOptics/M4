'''
@author: cs
'''

import os
import logging
import h5py
import numpy as np
from astropy.io import fits as pyfits
from m4.ground.interferometer_converter import InterferometerConverter
from m4.analyzer_iffunctions import AnalyzerIFF

class Noise():
    '''
    '''

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('NOISE:')
        self._ic = InterferometerConverter()



    def noise_acquisition_and_analysis(self, n_frame):
        file_path = self._acquisition(n_frame)
        rms_vector = self._analysis(file_path)
        return rms_vector

    def _acquisition(self, n_frame):
        file_path = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions/hdf5'
        return file_path

    def _createAndSaveCubeFromH5Data(self, file_path):
        self._cubeNoise = None
        for i in range(0, 4000):
            name = 'img_%04d.h5' %i
            file_name = os.path.join(file_path, name)
            image = self._ic.from4D(file_name)
            if self._cubeNoise is None:
                self._cubeNoise = image
            else:
                self._cubeNoise = np.ma.dstack((self._cubeNoise, image))

        fits_file_name = os.path.join(file_path, 'misure.fits')
        pyfits.writeto(fits_file_name, self._cubeNoise.data)
        pyfits.append(fits_file_name, self._cubeNoise.mask.astype(int))
        return self._cubeNoise

    def _readCube(self, file_path):
        fits_file_name = os.path.join(file_path, 'misure.fits')
        hduList = pyfits.open(fits_file_name)
        self._cubeNoise = np.ma.masked_array(hduList[0].data,
                                  hduList[1].data.astype(bool))
        return self._cubeNoise

    def _analysis(self):
        an = AnalyzerIFF.loadTestMeasureFromFits('hdf5') 

        pass
            