'''
@author: cs
'''

import os
import copy
import logging
import h5py
import numpy as np
from astropy.io import fits as pyfits
from m4.ground import tracking_number_folder
from m4.ground.configuration import Configuration
from m4.utils.img_redux import TipTiltDetrend
from m4.ground.interferometer_converter import InterferometerConverter
from m4.influence_functions_maker import IFFunctionsMaker
from m4.analyzer_iffunctions import AnalyzerIFF
from m4.utils.zernike_on_m_4 import ZernikeOnM4


class Noise():
    '''
    Classe per la valutazione del rumore

    HOW TO USE IT:
    from m4.noise_functions import Noise
    n = Noise()
    #acquisizione dati e analisi dalla cartella hdf5
    tt = n.noise_analysis_from_hdf5_folder(tidy_or_shuffle, template) 
    #analisi di più cartelle di dati
    rms_medio, tip, tilt, n_temp = n.different_template_analyzer(tt_list)
    '''

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('NOISE:')
        self._ic = InterferometerConverter()
        self._zOnM4 = ZernikeOnM4()
        self._ttd = TipTiltDetrend()
        self._numberOfNoiseIma = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return os.path.join(Configuration.OPD_DATA_FOLDER,
                            "Noise")

    ### IFF ##
    def noise_acquisition_and_analysis(self, numberOfNoiseIma, template=None):
        '''
        Funzione per la simulazione dell'acquisizione dati: una cartella per
        ogni template.

        args:
            template = np.array composed by 1 and -1
            numberOfNoiseIma = numero di immagini di rumore
                                (dalla cartella Noise/hdf5) da usare
                                NOTA: dipende dal sandbox.provaAcquisitionNoise ed è pari a
                                        n_temp * n_pushPull * n_modes

        returns:
            rms_mean = rms mediato sul numero di modi usato nell'acquisione delle iff
        '''
        if template is None:
            template = np.array([1, -1, 1])
        else:
            template = template
        self._numberOfNoiseIma = numberOfNoiseIma
        destination_file_path = self._acquisitionWithIFF(template)
        self._cubeFromAnalysis = self._analysisFromIFF(destination_file_path)
        rms_mean, quad_tt = self.rmsFromCube(self._cubeFromAnalysis)
        self._saveResults(rms_mean, quad_tt, destination_file_path)
        return destination_file_path

    def _acquisitionWithIFF(self, template):
        '''
        args:
            template = np.array composed by 1 and -1

        returns:
            destination_file_path = location of data
        '''
        # Aggiungere un'acquisizione consecutiva di tot misure in una sola cartella
        from m4 import sandbox
        destination_file_path = sandbox.provaAcquisitionNoise(template)
        return destination_file_path

    def _analysisFromIFF(self, destination_file_path):
        '''
        args:
            destination_file_path = dove prendere le info per far partire l'analisi

        returns:
            _cubeFromAnalysis = cubo ottenuto dopo l'analisi delle iff
        '''
        an = AnalyzerIFF.loadTestMeasureFromFits(destination_file_path)
        self._cubeFromAnalysis = an.createCube()

        fits_file_name = os.path.join(destination_file_path, 'Cube.fits')
        pyfits.writeto(fits_file_name, self._cubeFromAnalysis.data)
        pyfits.append(fits_file_name, self._cubeFromAnalysis.mask.astype(int))
        return self._cubeFromAnalysis

    def _createAndSaveCubeFromH5Data(self, data_file_path, destination_file_path, device):
        ''' Funzione che usa il procedimento di acquisizione dati delle iff
        args:
            data_file_path = cartella delle immagini di rumore
            destination_file_path = dove vengono messe le iff 
            device = quale segmento

        return:
            _cubeNoise = cubo di immagini delle iff
        '''
        IF = IFFunctionsMaker(device)
        self._cubeNoise = None
        for i in range(0, self._numberOfNoiseIma):
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
    ### Fine IFF ###

    def _defAnalyzer(self, tidy_or_shuffle, template, actsVector=None, n_push_pull=None):
        '''
        arg:
            tidy_or_shuffle = (int) 0 per tidy, 1 per shuffle
            template = np.array composed by 1 and -1
            actsVector = vector of actuators or modes
            n_push_pull = (int) number of push pull

        returns:
            an = analyzer object
        '''
        an = AnalyzerIFF()
        if n_push_pull is None:
            an._nPushPull = 3
        else:
            an._nPushPull = n_push_pull
        an._template = template
        if actsVector is None:
            an._actsVector = np.arange(25)
            an._modeVector = np.copy(an._actsVector)
        else:
            an._actsVector = actsVector
        an._cmdAmplitude = np.ones(an._actsVector.shape[0])

        indexingList = []
        if tidy_or_shuffle == 0:
            for i in range(an._nPushPull):
                indexingList.append(an._modeVector)
        elif tidy_or_shuffle == 1:
            for j in range(an._nPushPull):
                np.random.shuffle(an._modeVector)
                indexingList.append(an._modeVector)
        an._indexingList = np.array(indexingList)
        return an


    def noise_analysis_from_hdf5_folder(self, tidy_or_shuffle, template, actsVector=None, n_push_pull=None):
        '''
        arg:
            tidy_or_shuffle = (int) 0 per tidy, 1 per shuffle
            template = np.array composed by 1 and -1
            actsVector = vector of actuators or modes
            n_push_pull = (int) number of push pull

        returns:
            tt = tracking number of measurement
        '''
        data_file_path = os.path.join(Noise._storageFolder(), 'hdf5')

        an = self._defAnalyzer(tidy_or_shuffle, template, actsVector, n_push_pull)

        store_in_folder = self._storageFolder()
        save = tracking_number_folder.TtFolder(store_in_folder)
        dove, tt = save._createFolderToStoreMeasurements()

        self._cubeFromAnalysis = an.createCubeFromImageFolder(data_file_path)
        fits_file_name = os.path.join(dove, 'Cube.fits')
        pyfits.writeto(fits_file_name, self._cubeFromAnalysis.data)
        pyfits.append(fits_file_name, self._cubeFromAnalysis.mask.astype(int))
        fits_file_name = os.path.join(dove, 'template.fits')
        pyfits.writeto(fits_file_name, an._template)

        rms_mean, quad_mean = self.rmsFromCube(self._cubeFromAnalysis)
        self._saveResults(rms_mean, quad_mean, dove)
        return tt

    def rmsFromCube(self, cube_to_process):
        '''
        args:
            cube_to_process = cube generated by the analyzer_iffunctions

        returns:
            rms_mean = rms mediato sul numero di modi usato nell'acquisione delle iff
            tip = tip mediato sul numero di modi usato nell'acquisione delle iff
            tilt = tilt mediato sul numero di modi usato nell'acquisione delle iff
        '''
        rms_list = []
        coef_tilt_list = []
        coef_tip_list = []
        quad_list = []
        for i in range(cube_to_process.shape[2]):
            image = self._imageExtender(cube_to_process[:,:,i])
            coef, mat = self._zOnM4.zernikeFit(image,
                                               np.array([2, 3]))
            image_ttr = self._ttd.ttRemoverFromCoeff(coef, image)
            rms = image_ttr.std()
            rms_list.append(rms)
            coef_tip_list.append(coef[0])
            coef_tilt_list.append(coef[1])
            quad = np.sqrt(coef[0]**2 + coef[1]**2)
            quad_list.append(quad)
        rms_vector = np.array(rms_list)
        tip = np.array(coef_tip_list).mean()
        tilt = np.array(coef_tilt_list).mean()
        quad_tt = np.array(quad_list).mean()
        rms_mean = np.mean(rms_vector)
        return rms_mean, quad_tt

    def _imageExtender(self, cube_element):
        ''' Funzione usata per estendere le dimenzioni delle immagini acquisite a
        quelle della pupilla su cui costruisco gli zernike
        args:
            cube_element = an image from cube image

        returns:
            image = cube_element esteso alle dimensioni della pupilla degli zernike
        '''
        dim_y = (2 * self._zOnM4._zg.getRadius() - cube_element.shape[0]).astype(int) #512-500
        vv = np.ma.masked_array(np.zeros((dim_y, cube_element.shape[1])),
                                mask=np.ones((dim_y, cube_element.shape[1])).astype(bool)) #496
        dim_x = (2 * self._zOnM4._zg.getRadius() - cube_element.shape[1]).astype(int)   #512-496
        vv2 = np.ma.masked_array(np.zeros(((2 * self._zOnM4._zg.getRadius()).astype(int), dim_x)),
                                 mask=np.ones(((2 * self._zOnM4._zg.getRadius()).astype(int), dim_x)).astype(bool))
        pp = np.ma.append(cube_element, vv, axis=0)
        image = np.ma.append(pp, vv2, axis=1)
        return image

    def _saveResults(self, rms_mean, quad_mean, destination_file_path):
        ''' Save results as text file
        '''
        fits_file_name = os.path.join(destination_file_path, 'results.txt')
        file = open(fits_file_name, 'w+')
        file.write('%e %e' %(rms_mean, quad_mean))
        file.close()

    def _readResultsFromTxt(self, tt):
        ''' Read results from text file and return a string
        '''
        store_in_folder = Noise._storageFolder()
        file_path = os.path.join(store_in_folder, tt)
        file_name = os.path.join(file_path, 'results.txt')
        file = open(file_name, 'r')
        contents = file.read() #restituisce la riga come stringa
        return contents

    def _readCube(self, tt):
        '''
        args:
            tt = tracking number of measurement

        return:
            _cubeFromAnalysis = cubo ottenuto dopo l'analisi delle iff
        '''
        store_in_folder = Noise._storageFolder()
        file_path = os.path.join(store_in_folder, tt)
        fits_file_name = os.path.join(file_path, 'Cube.fits')
        hduList = pyfits.open(fits_file_name)
        self._cubeFromAnalysis = np.ma.masked_array(hduList[0].data,
                                  hduList[1].data.astype(bool))
        return self._cubeFromAnalysis

    def _readTempFromInfoFile(self, tt):
        '''
        args:
            tt = tracking number of measurement

        return:
            n_temp = (int) lunghezza del vettore template
        '''
        store_in_folder = Noise._storageFolder()
        file_path = os.path.join(store_in_folder, tt)
        fits_file_name = os.path.join(file_path, 'Info.fits')
        hduList= pyfits.open(fits_file_name)
        n_temp = hduList[1].data.shape[0]
        return n_temp

    def _readTempFromFits(self, tt):
        '''
        args:
            tt = tracking number of measurement

        returns:
            n_temp = (int) lunghezza del vettore template
        '''
        store_in_folder = Noise._storageFolder()
        file_path = os.path.join(store_in_folder, tt)
        fits_file_name = os.path.join(file_path, 'template.fits')
        hduList= pyfits.open(fits_file_name)
        n_temp = hduList[0].data.shape[0]
        return n_temp

    def different_template_analyzer(self, tt_list):
        '''
        args:
            tt_list = list of tracking number to analyze

        return:
            rms_medio = vettore degli rms medi (uno per ogni cartella di dati)
            n_tempo = vettore delle lunghezza dei template usati
        '''
        rms_list = []
        quad_list = []
        n_temp_list = []
        for tt in tt_list:
            cube = self._readCube(tt)
            n_temp = self._readTempFromFits(tt)
            #n_temp = self._readTempFromInfoFile(tt)
            rms, quad = self.rmsFromCube(cube)
            rms_list.append(rms)
            quad_list.append(quad)
            n_temp_list.append(n_temp)

        rms_medio = np.array(rms_list)
        quad = np.array(quad_list)
        n_temp = np.array(n_temp_list)
        return rms_medio, quad, n_temp
    ### tt_list ###
    # measurementFolder ='/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Noise'
    # list= os.listdir(measurementFolder); list.sort()
    # tt_list = list[1:len(list)-2]
    ### PLOT ###
    # plot(n_temp, rms, '-o'); plt.xlabel('n_temp'); plt.ylabel('rms_medio')

    ### FUNZIONE DI STRUTTURA ###
    def analysis_whit_structure_function(self, tau_vector):
        '''
        args:
            tau_vector = vector of tau to use
            4000 = total number of image in hdf5

        returns:
            rms_medio =
        '''
        data_file_path = os.path.join(Noise._storageFolder(), 'hdf5')
        list = os.listdir(data_file_path)
        image_number = len(list) -1
        #tau_vector = np.arange(80)+1
        i_max = np.int((image_number - tau_vector[tau_vector.shape[0]-1]) /
                       (tau_vector[tau_vector.shape[0]-1] * 2))
        if i_max <= 20:
            raise OSError('tau = %s too large' %tau_vector[tau_vector.shape[0]-1])
        rms_medio_list = []
        quad_med_list = []
        for j in range(tau_vector.shape[0]):
            dist = tau_vector[j]
            rms_list = []
            quad_list = []
            for i in range(i_max):
                k = i * dist * 2
                name = 'img_%04d.h5' %k
                file_name = os.path.join(data_file_path, name)
                image_k = self._ic.from4D(file_name)
                name = 'img_%04d.h5' %(k+dist)
                file_name = os.path.join(data_file_path, name)
                image_dist = self._ic.from4D(file_name)

                image_diff = image_k - image_dist
                image = self._imageExtender(image_diff)
                zernike_coeff_array, mat = self._zOnM4.zernikeFit(image,
                                                                  np.array([2, 3]))
                image_ttr = self._ttd.ttRemoverFromCoeff(zernike_coeff_array, image)
                quad = np.sqrt(zernike_coeff_array[0]**2 + zernike_coeff_array[1]**2)

                rms = image_ttr.std()
                rms_list.append(rms)
                quad_list.append(quad)
            rms_vector = np.array(rms_list)
            aa = rms_vector.mean()
            rms_medio_list.append(aa)
            quad_med_list.append(np.array(quad_list).mean())
        rms_medio = np.array(rms_medio_list)
        quad_med = np.array(quad_med_list)

        return rms_medio, quad_med
    # plot(tau_vector, rms, '-o'); plt.xlabel('tau'); plt.ylabel('rms_medio')
    