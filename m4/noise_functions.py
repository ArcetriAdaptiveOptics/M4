'''
@author: cs
'''

import os
import logging
import h5py
import numpy as np
from astropy.io import fits as pyfits
from m4.ground.interferometer_converter import InterferometerConverter
from m4.influence_functions_maker import IFFunctionsMaker
from m4.analyzer_iffunctions import AnalyzerIFF

class Noise():
    '''
    '''

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('NOISE:')
        self._ic = InterferometerConverter()



    def noise_acquisition_and_analysis(self,):
        destination_file_path = self._acquisition()
        self._cubeFromAnalysis = self._analysis(destination_file_path)
        rms_vector = self.rmsFromCube(self._cubeFromAnalysis)
        return rms_vector

    def _acquisition(self):
        from m4 import sandbox
        destination_file_path = sandbox.provaAcquisitionNoise()
        return destination_file_path

    def rmsFromCube(self, cube_to_process):
        rms_list = []
        for i in range(cube_to_process.shape[2]):
            rms = cube_to_process[:,:,i].std()
            rms_list.append(rms)
        rms_vector = np.array(rms_list)
        return rms_vector

    def _createAndSaveCubeFromH5Data(self, data_file_path, destination_file_path, device):
        IF = IFFunctionsMaker(device)
        self._cubeNoise = None
        for i in range(0, 500):
            name = 'img_%04d.h5' %i
            file_name = os.path.join(data_file_path, name)
            image = self._ic.from4D(file_name)
            if self._cubeNoise is None:
                self._cubeNoise = image
            else:
                self._cubeNoise = np.ma.dstack((self._cubeNoise, image))

        who, tt_cmdH, acts_vector, cmd_matrix, \
            amplitude, n_push_pull, indexingList, template = IF.loadInfoFromFits(destination_file_path)

        fits_file_name = os.path.join(destination_file_path, 'misure.fits')
        header = pyfits.Header()
        header['NPUSHPUL'] = n_push_pull
        header['WHO'] = who
        header['TT_CMDH'] = tt_cmdH
        pyfits.writeto(fits_file_name, acts_vector, header)
        pyfits.append(fits_file_name, cmd_matrix, header)
        pyfits.append(fits_file_name, amplitude, header)
        pyfits.append(fits_file_name, indexingList, header)
        pyfits.append(fits_file_name, self._cubeNoise.data, header)
        pyfits.append(fits_file_name, self._cubeNoise.mask.astype(int), header)
        pyfits.append(fits_file_name, template, header)
        return self._cubeNoise

    def _readCube(self, file_path):
        fits_file_name = os.path.join(file_path, 'Cube.fits')
        hduList = pyfits.open(fits_file_name)
        self._cubeFromAnalysis = np.ma.masked_array(hduList[0].data,
                                  hduList[1].data.astype(bool))
        return self._cubeFromAnalysis

    def _analysis(self, destination_file_path):
        an = AnalyzerIFF.loadTestMeasureFromFits(destination_file_path)
        self._cubeFromAnalysis = an.createCube()

        fits_file_name = os.path.join(destination_file_path, 'Cube.fits')
        pyfits.writeto(fits_file_name, self._cubeFromAnalysis.data)
        pyfits.append(fits_file_name, self._cubeFromAnalysis.mask.astype(int))
        return self._cubeFromAnalysis
