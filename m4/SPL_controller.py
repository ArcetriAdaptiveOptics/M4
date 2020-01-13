'''
@author: cs
'''

import os
import logging
import glob
import time
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
        '''
        args:
            lambda_vector = vector of wavelengths (between 400/700 nm)
                            [np.array]
            acqu4d = 1 to perform the data acquisition process
                     0 to skip data acquisition
            an = 1 to perform the data analysis process
                 0 to skip data analysis

        returns:
                piston = piston value
        '''
        if (acq4d is None and an is None):
            self._exptime, self._acq4d, self._an = self._setParameters(0.7, 1, 1)
            tt = self.acquire(lambda_vector)
            piston = self.analyzer(tt)

        elif (acq4d==1 and an==0):
            self._exptime, self._acq4d, self._an = self._setParameters(0.7, 1, 0)
            tt = self.acquire(lambda_vector)

        elif (acq4d==0 and an==1):
            self._exptime, self._acq4d, self._an = self._setParameters(0.7, 0, 1)
            tt = input('Tracking number of measurement to analyze:')
            piston = self.analyzer(tt)

        else:
            raise OSError('Combination not permitted')

        return piston

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
        '''
        arg:
            lambda_vector = vector of wavelengths (between 400/700 nm)

        returns:
                tt = tracking number of measurements
        '''
        store_in_folder = self._storageFolder()
        save = tracking_number_folder.TtFolder(store_in_folder)
        self._dove, self._tt = save._createFolderToStoreMeasurements()
        fits_file_name = os.path.join(self._dove, 'lambda_vector.fits')
        pyfits.writeto(fits_file_name, lambda_vector)

        ## find PSF position ##
        self._filter.setWavelength(600)
        ExposureTimeAbs = 1.5 * self._exptime * 1e6
        self._camera.setExposureTime(ExposureTimeAbs)
        self._camera.Timeout = 30
        reference_image = self._camera.acquireFrame() #double snapshot
        #mettergli la maschera gpixmask
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

        self._logger.info('Acquisition of frames')
        for i in range(lambda_vector.shape[0]):
            self._filter.setWavelength(lambda_vector[i])
            self._camera.ExposureTime = self._exptime * expgain[i] * 1e6
            time.sleep(3 * ExposureTimeAbs / 1e6)
            image = self._camera.acquireFrame() #single snapshot
            #mettergli la maschera gpixmask
            crop = self._preProcessing(image)

            file_name = 'image_%dnm.fits' %lambda_vector[i]
            self._saveCameraFrame(file_name, crop)

        return self._tt

    def _baricenterCalculator(self, reference_image):
        ''' Calculate the peak position of the image
        args:
            reference_image = camera frame
        '''
        counts, bin_edges = np.histogram(reference_image)
        thr = 5 * bin_edges[np.where(counts == max(counts))]
        idx = np.where(reference_image < thr)
        reference_image[idx] = 0
        cy, cx = scin.measurements.center_of_mass(reference_image)
        self._baricenterCoord = [np.int(cy), np.int(cx)]


    def _preProcessing(self, image):
        ''' Cut the images around the pick position
        args:
            image = camera frame

        returns:
            crop = cut camera frame
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
        img = image - bkg #va mascherata con gpixmask
        crop = img[self._baricenterCoord[0]-ycrop:
                   self._baricenterCoord[0]+ycrop,
                   self._baricenterCoord[1]-xcrop:
                   self._baricenterCoord[1]+xcrop]
        return crop

    def _saveCameraFrame(self, file_name, frame):
        ''' Save camera frames in SPL/TrackingNumber
        '''
        fits_file_name = os.path.join(self._dove, file_name)
        pyfits.writeto(fits_file_name, frame.data)



    def analyzer(self, tt=None):
        ''' Analyze measurement data and compare it with synthetic data
        args:
            tt = tracking number of the measurement data

        returns:
            piston = piston value
        '''
        if tt is None:
            tt = self._tt
            dove = self._dove
        else:
            tt = tt
            dove = os.path.join(Configuration.LOG_ROOT_FOLDER, 'SPL', tt)

        self._logger.info('Analysis of tt = %s', tt)

        lambda_path = os.path.join(dove, 'lambda_vector.fits')
        hduList = pyfits.open(lambda_path)
        lambda_vector = hduList[0].data

        cube, cube_normalized = self._readCameraFrames(tt)
        img = np.sum(cube_normalized, 2)
        pick = self._newThr(img)
        matrix = np.zeros((pick[3]-pick[2] + 1, lambda_vector.shape[0])) # 150 pixel
        matrix_smooth = np.zeros((pick[3]-pick[2] + 2, lambda_vector.shape[0]))
        crop_frame_cube = None
        for i in range(lambda_vector.shape[0]):
            frame = cube[:,:,i]
            crop_frame = frame[pick[0]:pick[1], pick[2]:pick[3] + 1]
            if np.max(crop_frame) > 4000:
                print("**************** WARNING: saturation detected!")
            if crop_frame_cube is None:
                crop_frame_cube = crop_frame
            else:
                crop_frame_cube = np.dstack((crop_frame_cube, crop_frame))

            y = np.sum(crop_frame, 0)
            area = np.sum(y[:])
            y_norm = y / area
            if i == 0:
                mm = 1.2 * np.max(y_norm)
                matrix[:,i] = mm
            else:
                matrix[:,i] = y_norm

            w = self._smooth(y_norm, 4)
            matrix_smooth[:, i] = w

        piston = self.templateComparison(matrix, lambda_vector)

        return piston



    def _readCameraFrames(self, tt):
        ''' Read images in a specific tracking number and return the cube
        args:
            tt = tracking number

        returns:
            cube = [pixels, pixels, n_frames=lambda]
            cube_normalized = creating with normalized images,
                                [pixels, pixels, n_frames=lambda]
        '''
        dove = os.path.join(Configuration.LOG_ROOT_FOLDER, 'SPL', tt)

        path_list = []
        for f in  glob.iglob(os.path.join(dove, 'image_*.fits')): 
            path_list.append(f)

        path_list.sort()
        cube = None
        cube_normalized = None
        for i in range(len(path_list)):
            hduList = pyfits.open(path_list[i])
            image = hduList[0].data
            if cube is None:
                cube = image
            else:
                cube = np.dstack((cube, image))
            image_normalized = image / np.sum(image)
            if cube_normalized is None:
                cube_normalized = image_normalized
            else:
                cube_normalized = np.dstack((cube_normalized, image_normalized))
        return cube, cube_normalized

    def _newThr(self, img):
        ''' Calculate the peak position of the image '''
        counts, bin_edges = np.histogram(img)
        edges = (bin_edges[2:] + bin_edges[1:len(bin_edges)-1]) / 2
        thr = 5 * edges[np.where(counts == max(counts))]
        idx = np.where(img < thr)
        img[idx] = 0
        cy, cx = scin.measurements.center_of_mass(img)
        baricenterCoord = [np.int(round(cy)), np.int(round(cx))]
        pick = [baricenterCoord[0]-25, baricenterCoord[0]+25,
                baricenterCoord[1]-75, baricenterCoord[1]+75]
        return pick


### FUNZIONE PER SMOOTH ###

    def _smooth(self, a, WSZ):
        ''''
        args:
            a = NumPy 1-D array containing the data to be smoothed
            WSZ = smoothing window size needs, which must be odd number,
                 as in the original MATLAB implementation '''
        out0 = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ
        r = np.arange(1,WSZ-1,2)
        start = np.cumsum(a[:WSZ-1])[::2]/r
        stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
        return np.concatenate((  start , out0, stop  ))

### FINE ###


    def templateComparison(self, matrix, lambda_vector):
        ''' Compare the matrix obtained from the measurements with
            the one recreated with the synthetic data in tn_fringes
        args:
            matrix = [pixels, lambda]
            lambda_vector = vector of wavelengths

        returns:
            piston = piston value
        '''
        tn_fringes = '20181108_1'
        #tn_fringes = '20181026'
        self._logger.debug('Template Comparison with data in tn_fringes = %s',
                           tn_fringes)
        dove = os.path.join(Configuration.LOG_ROOT_FOLDER, 'SPL/Fringes',
                            tn_fringes)
        dove_delta = os.path.join(dove, 'Differential_piston.fits')
        hduList = pyfits.open(dove_delta)
        delta = hduList[0].data
        dove_lambda_synth = os.path.join(dove, 'Lambda.fits')
        hduList = pyfits.open(dove_lambda_synth)
        lambda_synth_from_data = hduList[0].data * 1e9
        lambda_synth = self._myLambaSynth(lambda_synth_from_data)

        idx = np.isin(lambda_synth, lambda_vector)

        Qm = matrix - np.mean(matrix)
        #creare la matrice di Giorgio della giusta dimenzione
        F = None
        for i in range(1, delta.shape[0]):
            file_name = os.path.join(dove, 'Fringe_%05d.fits' %i)
            hduList = pyfits.open(file_name)
            fringe = hduList[0].data
            fringe_selected = fringe[:,idx]
            if F is None:
                F = fringe_selected
            else:
                F = np.dstack((F, fringe_selected))
        Qt = F - np.mean(F)

        R = np.zeros(delta.shape[0]-1)
        for i in range(delta.shape[0]-1):
            R[i] = np.sum(np.sum(Qm[:,:]*Qt[:,:,i])) / ((np.sum(np.sum(Qm[:,:]**2)))**5
                                                    * (np.sum(np.sum(Qt[:,:,i]**2)))**5)

        idp = np.where(R==max(R))
        piston = delta[idp]

        return piston

    def _myLambaSynth(self, lambda_synth_from_data):
        ''' Transforms its values into integers
        '''
        my_lambda = np.zeros(lambda_synth_from_data.shape[0])
        for j in range(lambda_synth_from_data.shape[0]):
            bb = np.int(round(lambda_synth_from_data[j]))
            my_lambda[j] = bb
        return my_lambda
