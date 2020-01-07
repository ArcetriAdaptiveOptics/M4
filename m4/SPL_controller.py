'''
@author: cselmi
'''

import os
import logging
import glob
import numpy as np
import scipy.ndimage as scin
from m4.ground import tracking_number_folder
from astropy.io import fits as pyfits
from m4.ground.configuration import Configuration

class SPL():

    def __init__(self, filter_obj, camera_obj):
        """The constructor """
        self._logger = logging.getLogger('SPL_CONTROLLER:')
        self._filter = filter_obj
        self._camera = camera_obj

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "SPL")

    def measurement(self, lambda_vector, acq4d=None, an=None):
        ''' Deve contenere tutta la procedura di acquisizione ed analisi dei dati
        '''
        if (acq4d is None and an is None):
            self._exptime, self._acq4d, self._an = self._setParameters(0.7, 1, 1)
            self.acquire(lambda_vector)
            self.analyzer()

        elif (acq4d==1 and an==0):
            self._exptime, self._acq4d, self._an = self._setParameters(0.7, 1, 0)
            self.acquire(lambda_vector)

        else:
            raise OSError('Combination not permitted')

    def _setParameters(self, exptime, acq4d, an):
        '''
        args:
            exptime = exposure time
            acq4d = 0 to not acquire data, 1 to acquire data
            an = 0 to not analyze data, 1 for the analysis
        '''
        self._exptime = exptime
        self._acq4d = acq4d
        self._an = an
        return self._exptime, self._acq4d, self._an


# lambda_vector = np.arange(530,730,10)
    def acquire(self, lambda_vector):
        ''' Acquisizione dati
        arg:
            lambda_vector = vector of wavelengths (between 400/700 nm)
        '''
        store_in_folder = self._storageFolder()
        save = tracking_number_folder.TtFolder(store_in_folder)
        self._dove, self._tt = save._createFolderToStoreMeasurements()

        ## find PSF position ##
        self._filter.setWavelength(600)
        ExposureTimeAbs = 1.5 * self._exptime * 1e6
        self._camera.setExposureTime(ExposureTimeAbs)
        self._camera.Timeout = 30
        reference_image = self._camera.acquireFrame() #mettergli la maschera gpixmask
        if np.max(reference_image) > 4000:
            print("**************** WARNING: saturation detected!")
        # calcolo il baricentro
        self._baricenterCalculator(reference_image)


        expgain= np.ones(lambda_vector.shape[0])
        expgain[np.where(lambda_vector < 500)] = 8
        expgain[np.where(lambda_vector < 530)] = 4
        expgain[np.where(lambda_vector > 680)] = 3
        expgain[np.where(lambda_vector > 700)] = 8
        expgain[np.where(lambda_vector >= 720)] = 9


        crop_cube = None
        crop_cube_normalized = None
        for i in range(lambda_vector.shape[0]):
            self._filter.setWavelength(lambda_vector[i])
            self._camera.ExposureTime = self._exptime * expgain(i) * 1e6
            image = self._camera.acquireFrame()
            crop = self._preProcessing(image)
            crop_normalized = crop / np.sum(crop)

            file_name = 'image_%dnm.fits' %lambda_vector[i]
            self._saveCameraFrame(file_name, crop)

            if crop_cube is None:
                crop_cube = crop
            else:
                crop_cube = np.dstack((crop_cube, crop))
            if crop_cube_normalized is None:
                crop_cube_normalized = crop_normalized
            else:
                crop_cube_normalized = np.dstack((crop_cube_normalized, crop_normalized))

        return crop_cube, crop_cube_normalized

    def _baricenterCalculator(self, reference_image):
        ''' Calcola le coordinate del baricentro dell'immagine
        '''
        counts, bin_edges = np.histogram(reference_image)
        thr = 5 * bin_edges[np.where(counts == max(counts))]
        idx = np.where(reference_image < thr)
        reference_image[idx] = 0
        cy, cx = scin.measurements.center_of_mass(reference_image)
        self._baricenterCoord = [np.int(cy), np.int(cx)]


    def _preProcessing(self, image):
        ''' Prepara i dati acquisiti per essere analizzati
        '''
        xcrop = 145 #150
        ycrop = 95 #100

        tmp = np.zeros((image.shape[0], image.shape[1]))
        tmp[self._baricenterCoord[0]-ycrop:
            self._baricenterCoord[0]+ycrop,
            self._baricenterCoord[1]-xcrop:
            self._baricenterCoord[1]+xcrop] = 1
        id_bkg = np.where(tmp == 0)
        bkg = np.ma.mean(image[id_bkg])
        img = image - bkg #va mascherata
        crop = img[self._baricenterCoord[0]-ycrop:
                   self._baricenterCoord[0]+ycrop,
                   self._baricenterCoord[1]-xcrop:
                   self._baricenterCoord[1]+xcrop]
        return crop

    def _saveCameraFrame(self, file_name, frame):
        ''' Salva le immagini acquisite dalla telecamera
        '''
        fits_file_name = os.path.join(self._dove, file_name)
        pyfits.writeto(fits_file_name, frame.data)



    def analyzer(self):
        ''' Analizza i dati e li confronta con quelli sintetici
        '''

    def _readCameraFrames(self, tt=None):
        ''' Legge le immagini in un determinato tracking number e ne restituisce il cubo
        '''
        if tt is None:
            tt = self._tt
            dove = self._dove
        else:
            tt = tt
            dove = os.path.join(Configuration.LOG_ROOT_FOLDER, 'SPL', tt)

        path_list = []
        for f in  glob.iglob(os.path.join(dove, 'image_*.fits')): 
            path_list.append(f)

        path_list.sort()
        cube = None
        for i in range(len(path_list)):
            hduList = pyfits.open(path_list[i])
            image = hduList[0].data
            if cube is None:
                cube = image
            else:
                cube = np.dstack((cube, image))
        return cube


