'''
@author: cs
'''

import os
import logging
import h5py
import numpy as np
from astropy.io import fits as pyfits
from m4.ground.configuration import Configuration
from m4.ground.interferometer_converter import InterferometerConverter
from m4.influence_functions_maker import IFFunctionsMaker
from m4.analyzer_iffunctions import AnalyzerIFF
from m4.utils.zernike_on_m_4 import ZernikeOnM4


class Noise():
    '''
    '''

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('NOISE:')
        self._ic = InterferometerConverter()
        self._zOnM4 = ZernikeOnM4()

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return os.path.join(Configuration.OPD_DATA_FOLDER,
                            "Noise")

    def noise_acquisition_and_analysis(self,):
        destination_file_path = self._acquisition()
        self._cubeFromAnalysis = self._analysis(destination_file_path)
        rms_mean, tip, tilt = self.rmsFromCube(self._cubeFromAnalysis)
        return rms_mean

    def _acquisition(self):
        from m4 import sandbox
        destination_file_path = sandbox.provaAcquisitionNoise()
        return destination_file_path

    def rmsFromCube(self, cube_to_process):
        rms_list = []
        coef_tilt_list = []
        coef_tip_list = []
        for i in range(cube_to_process.shape[2]):
            rms = cube_to_process[:,:,i].std()
            image = self._imageExtender(cube_to_process[:,:,i])
            coef, mat = self._zOnM4.zernikeFit(image,
                                               np.array([2, 3])) #non tornano i punti
            rms_list.append(rms)
            coef_tip_list.append(coef[0])
            coef_tilt_list.append(coef[1])
        rms_vector = np.array(rms_list)
        tip = np.array(coef_tip_list).mean()
        tilt = np.array(coef_tilt_list).mean()
        rms_mean = np.mean(rms_vector)
        return rms_mean, tip, tilt

    def _imageExtender(self, cube_element):
        vv = np.ma.masked_array(np.zeros((12, 496)), mask=np.ones((12, 496)).astype(bool))
        vv2 = np.ma.masked_array(np.zeros((512, 16)), mask=np.ones((512, 16)).astype(bool))
        pp = np.ma.append(cube_element, vv, axis=0)
        image = np.ma.append(pp, vv2, axis=1)
        return image

    def _createAndSaveCubeFromH5Data(self, data_file_path, destination_file_path, device):
        IF = IFFunctionsMaker(device)
        self._cubeNoise = None
        for i in range(0, 1876):
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

    def _readCube(self, tt):
        store_in_folder = Noise._storageFolder()
        file_path = os.path.join(store_in_folder, tt)
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

    def _readTempFromInfoFile(self, tt):
        store_in_folder = Noise._storageFolder()
        file_path = os.path.join(store_in_folder, tt)
        fits_file_name = os.path.join(file_path, 'Info.fits')
        hduList= pyfits.open(fits_file_name)
        n_temp = hduList[1].data.shape[0]
        return n_temp

    def different_template_analyzer(self, tt_list):
        rms_list = []
        n_temp_list = []
        for tt in tt_list:
            cube = self._readCube(tt)
            n_temp = self._readTempFromInfoFile(tt)
            rms, tip, tilt = self.rmsFromCube(cube)
            rms_list.append(rms)
            n_temp_list.append(n_temp)

        rms_medio = np.array(rms_list)
        n_temp = np.array(n_temp_list)
        return rms_medio, n_temp
    ### tt_list ###
    # measurementFolder ='/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Noise'
    # list= os.listdir(measurementFolder); list.sort()
    # tt_list = list[1:len(list)-1]
    ### PLOT ###
    # plot(n_temp, rms, '-o'); plt.xlabel('n_temp'); plt.ylabel('rms_medio')
