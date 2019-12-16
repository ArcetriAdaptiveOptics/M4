'''
@author: cs
'''

import os
import logging
import h5py
from astropy.io import fits as pyfits
import numpy as np
from m4.ground.interferometer_converter import InterferometerConverter
from m4.influence_functions_maker import IFFunctionsMaker
from m4.utils.roi import ROI
from m4.utils.img_redux import TipTiltDetrend
from m4.ground.configuration import Configuration


class AnalyzerIFF():
    '''
    This class analyzes the measurements made through the IFF class by
    generating the cube of measurements and calculating interaction matrix
    and reconstructor.

    HOW TO USE IT:
    from m4.analyzer_iffunctions import AnalyzerIFF
    fileName = os.path.join(".../IFFunctions", tt)
    an = AnalyzerIFF.loadInfoFromTtFolder(fileName)
    cube = an.createCube(tiptiltDetrend = None, phaseAmbiguity = None)
    an.saveCubeAsFits(cubeName)
    '''

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('IFF_ANALYZER:')
        self._ic = InterferometerConverter()
        self._indexingList = None
        self._cube = None
        self._rec = None
        self._intMat = None
        self._analysisMask = None
        self._cubeMeasure = None

    def _ttData(self):
        """ Allow to obtain the tracking number from the path """
        split = os.path.split(self._h5Folder)
        self._tt = split[1]
        return self._tt

    @staticmethod
    def loadInfoFromTtFolder(h5Folder):
        """ Creates the object using information about path measurements

            Args:
                h5Folder = measurement folder path

            Return:
                theObject = analyzerIFF class object
        """
        theObject = AnalyzerIFF()
        theObject._h5Folder = h5Folder
        a = IFFunctionsMaker.loadInfoFromFits(h5Folder)
        theObject._who = a[0]
        theObject._tt = a[1]
        theObject._actsVector = a[2]
        theObject._cmdMatrix = a[3]
        theObject._cmdAmplitude = a[4]
        theObject._nPushPull = a[5]
        theObject._indexingList = a[6]
        return theObject



    def getCube(self):
        """ Returns the cube of measurements """
        return self._cube

    def getIFShape(self):
        return self.getCube()[:,:,0].shape

    def getMasterMask(self):
        """ Returns the mask of the cube """
        aa = np.sum(self._cube.mask.astype(int), axis=2)
        master_mask = np.zeros(aa.shape, dtype=np.bool)
        master_mask[np.where(aa > 0)] = True
        return master_mask


    def _indexReorganization(self):
        """ Returns the index position """
        indv = np.array(self._indexingList)
        where = []
        for ind in self._actsVector:
            for j in range(self._nPushPull):
                a = np.where(indv[j] == ind)
                where.append(a)
        return where

    def _amplitudeReorganization(self, indexing_input, indexing_list,
                                 amplitude, n_push_pull):
        '''
            Args:
                indexing_input = vector of selected modes for carrying out
                                 the influence functions

                indexing_list = tuple indicating how the modes were applied

                amplitude = amplitude of applied modes

            Returns:
                vect = vector (amp.shape x n_push_pull.shape) with the
                        amplitudes ordered in the same way as indexing_list

        '''
        where = []
        for i in indexing_input:
            for j in range(n_push_pull):
                a = np.where(indexing_list[j] == i)
                where.append(a)
        where = np.array(where)
        vect = np.zeros(amplitude.shape[0]*n_push_pull)

        for i in range(amplitude.shape[0]):
            for k in range(n_push_pull):
                p = n_push_pull * i + k
                indvect = where[p][0][0]+ indexing_input.shape[0] * k
                vect[indvect] = amplitude[i]
        return vect

    def createTestCubeMeasure_forTestDataInOpdImages(self):
        tt = '20191210_143625'
        fold = os.path.join(Configuration.CALIBRATION_ROOT_FOLDER, 'OPDImages', tt, 'hdf5')

        hduList = pyfits.open(os.path.join(fold, 'ampVect.fits'))
        self._cmdAmplitude = hduList[0].data
        hduList = pyfits.open(os.path.join(fold, 'modeRange.fits'))
        self._actsVector = hduList[0].data
        self._nPushPull = 1 #per runa sarebbe il vettore array([ 1, -1,  1]), ossia 1.5 per me!
        indList = []
        for i in range(self._nPushPull):
            indList.append(self._actsVector)
        self._indexingList = np.array(indList)
        
        self._h5Folder = fold
        self._who = 'Bho'
#         for i in range(2675):
#             file_name = os.path.join(fold, 'img_%04d.h5') %i
#             ima = self._ic.from4D(file_name)
#             if self._cubeMeasure is None:
#                 self._cubeMeasure = ima
#             else:
#                 self._cubeMeasure = np.ma.dstack((self._cubeMeasure, ima))
# 
#         cube_name= 'CubeMeasure.fits'
#         fits_file_name = os.path.join(fold, cube_name)
#         hduList = pyfits.open(fits_file_name)
#         cube_measure = np.ma.masked_array(hduList[0].data,
#                                   hduList[1].data.astype(bool))
#         self._cubeMeasure = cube_measure
        return self._cubeMeasure

        

    def createCube(self, tiptilt_detrend=None, phase_ambiguity=None):
        '''
            Args:
                ttDetrend = in the creation of the cube the images are reduced
                            removing tip tilt on the central segment

                phaseSolve=

            Returns:
                self._cube = cube
        '''
        cube_all_act = None
        self._ttData()
        self._logCubeCreation(tiptilt_detrend, phase_ambiguity)
        where = self._indexReorganization()
        misure_pos, misure_neg = self._splitMeasureFromFits(self._cubeMeasure)
        ampl_reorg = self._amplitudeReorganization(self._actsVector,
                                                   self._indexingList,
                                                   self._cmdAmplitude,
                                                   self._nPushPull)


        for i in range(self._actsVector.shape[0]):
            for k in range(self._nPushPull):
                p = self._nPushPull * i + k
                #where= self._indexReorganization()
                n = where[p][0][0]
                mis = k* self._indexingList.shape[1] + n

                img_pos = misure_pos[:,:,mis]
                img_neg = misure_neg[:,:,mis]
                img_if = (img_pos - img_neg) / (2 * ampl_reorg[mis])
                if tiptilt_detrend is None:
                    img_if = img_if
                else:
                    r = ROI()
                    roi = r.roiGenerator(img_if)
                    tt = TipTiltDetrend()
                    img_if = tt.tipTiltRemover(img_if, roi, 3)

                if_push_pull_kth = img_if-np.ma.median(img_if)

                if k == 0:
                    all_push_pull_act_jth = if_push_pull_kth
                else:
                    all_push_pull_act_jth = np.ma.dstack((all_push_pull_act_jth,
                                                          if_push_pull_kth))

            if self._nPushPull == 1:
                if_act_jth = all_push_pull_act_jth
            else:
                if_act_jth = np.ma.mean(all_push_pull_act_jth, axis=2)

            if cube_all_act is None:
                cube_all_act = if_act_jth
            else:
                cube_all_act = np.ma.dstack((cube_all_act, if_act_jth))
        self._cube = cube_all_act

        return self._cube
    
    def createCube_forTestDataInOpdImages(self):
        tt = '20191210_143625'
        fold = os.path.join(Configuration.CALIBRATION_ROOT_FOLDER, 'OPDImages', tt, 'hdf5')

        hduList = pyfits.open(os.path.join(fold, 'ampVect.fits'))
        self._cmdAmplitude = hduList[0].data
        hduList = pyfits.open(os.path.join(fold, 'modeRange.fits'))
        self._actsVector = hduList[0].data
        hduList = pyfits.open(os.path.join(fold, 'template.fits'))
        vector_of_push_pull = hduList[0].data

        cube = None
        for i in range(self._actsVector.shape[0]):
            k = i * vector_of_push_pull.shape[0]
            imaList = []
            for j in range(vector_of_push_pull.shape[0]):
                l = k + j
                file_name = os.path.join(fold, 'img_%04d.h5') %l
                ima = self._ic.from4D(file_name)
                imaList.append(ima)

            image = imaList[0] + imaList[1] - imaList[2]

            if cube is None:
                cube = image 
            else:
                cube = np.ma.dstack((cube, image))

                
        return cube


    def _splitMeasureFromFits(self, misure):
        """ Divides the cube of measurements into two cubes

            Args:
                misure= cube of measurements

            Returns:
                misure_pos = positive measurements
                misure_neg = negative measurements
        """
        misure_pos = None
        misure_neg = None
        for j in range(misure.shape[2]):
            if j%2 == 0:
                if misure_pos is None:
                    misure_pos = misure[:,:,j]
                else:
                    misure_pos = np.ma.dstack((misure_pos, misure[:,:,j]))
            else:
                if misure_neg is None:
                    misure_neg = misure[:,:,j]
                else:
                    misure_neg = np.ma.dstack((misure_neg, misure[:,:,j]))

        return misure_pos, misure_neg

    def _logCubeCreation(self, tiptilt_detrend=None, phase_ambiguity=None):
        """ Use logging function to keep a record of
        how the cube was created """
        if (tiptilt_detrend is None and phase_ambiguity is None):
            self._logger.info('Creation of the IFF cube for %s. Location: %s',
                              self._who, self._tt)
            self._logger.debug('(Tip and tilt= ignored, \
                                    Phase ambiguity= ignored)')
        elif tiptilt_detrend is None:
            self._logger.info('Creation of the IFF cube for %s. Location: %s',
                              self._who, self._tt)
            self._logger.debug('(Tip and tilt= ignored, \
                                    Phase ambiguity= resolved)')
        elif phase_ambiguity is None:
            self._logger.info('Creation of the IFF cube for %s. Location: %s',
                              self._who, self._tt)
            self._logger.debug('(Tip and tilt= removed, \
                                     Phase ambiguity= ignored)')
        else:
            self._logger.info('Creation of the IFF cube for %s. Location: %s',
                              self._who, self._tt)
            self._logger.debug('(Tip and tilt= removed, \
                                    Phase ambiguity= resolved)')


    def saveCubeAsFits(self, cube_name):
        """
            Args:
                cube_name = name to save the cube
                            example 'Cube.fits'
        """
        tt = self._ttData()
        dove = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, tt)
        fits_file_name = os.path.join(dove, cube_name)
        header = pyfits.Header()
        header['NPUSHPUL'] = self._nPushPull
        header['WHO'] = self._who
        pyfits.writeto(fits_file_name, self._cube.data, header)
        pyfits.append(fits_file_name, self._cube.mask.astype(int))
        pyfits.append(fits_file_name, self._cmdAmplitude)
        pyfits.append(fits_file_name, self._actsVector)

    def saveCubeAsH5(self, cube_name):
        """
            Args:
                cube_name = name to save the cube
                            example 'Cube.h5'
        """
        tt = self._ttData()
        dove = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, tt)
        fits_file_name = os.path.join(dove, cube_name)
        hf = h5py.File(fits_file_name, 'w')
        hf.create_dataset('dataset_1', data=self._cube.data)
        hf.create_dataset('dataset_2', data=self._cube.mask.astype(int))
        hf.create_dataset('dataset_3', data=self._cmdAmplitude)
        hf.create_dataset('dataset_4', data=self._actsVector)
        hf.attrs['NPUSHPUL'] = self._nPushPull
        hf.attrs['WHO'] = self._who
        hf.close()


    @staticmethod
    def loadCubeFromFits(fits_file_name):
        """ Creates the object using information contained in Cube

            Args:
                fits_file_name = cube path

            Return:
                theObject = analyzerIFF class object
        """
        header = pyfits.getheader(fits_file_name)
        hduList = pyfits.open(fits_file_name)
        cube = np.ma.masked_array(hduList[0].data,
                                  hduList[1].data.astype(bool))
        acts_vector = hduList[3].data
        cmd_amplitude = hduList[2].data
        who = header['WHO']
        try:
            n_push_pull = header['NPUSHPUL']
        except KeyError:
            n_push_pull = 1

        theObject = AnalyzerIFF()
        theObject._actsVector = acts_vector
        theObject._cmdAmplitude = cmd_amplitude
        theObject._nPushPull = n_push_pull
        theObject._who = who
        theObject._cube = cube
        return theObject

    @staticmethod
    def loadCubeFromH5(file_name):
        """ Creates the object using information contained in Cube

            Args:
                file_name = cube path

            Return:
                theObject = analyzerIFF class object
        """
        theObject = AnalyzerIFF()
        hf = h5py.File(file_name, 'r')
        hf.keys()
        data1 = hf.get('dataset_1')
        data2 = hf.get('dataset_2')
        data3 = hf.get('dataset_3')
        data4 = hf.get('dataset_4')
        theObject._cube = np.ma.masked_array(np.array(data1),
                                             np.array(data2.astype(bool)))
        theObject._actsVector = np.array(data4)
        theObject._cmdAmplitude = np.array(data3)
        theObject._nPushPull = hf.attrs['NPUSHPUL']
        theObject._who = hf.attrs['WHO']
        hf.close()
        return theObject


    @staticmethod
    def loadTestMeasureFromFits(tt):
        """
        Creates the AnalyzerIFF object using specific informations
        """
        dove = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, tt)
        fits_file_name = os.path.join(dove, 'misure.fits')
        header = pyfits.getheader(fits_file_name)
        hduList = pyfits.open(fits_file_name)
        theObject = AnalyzerIFF()
        theObject._h5Folder = dove
        theObject._cubeMeasure = \
            np.ma.masked_array(hduList[4].data, hduList[5].data.astype(bool))
        theObject._actsVector = hduList[0].data
        theObject._cmdMatrix = hduList[1].data
        theObject._cmdAmplitude = hduList[2].data
        theObject._indexingList = hduList[3].data
        who = header['WHO']
        tt_cmdH = header['TT_CMDH']
        try:
            n_push_pull = header['NPUSHPUL']
        except KeyError:
            n_push_pull = 1
        theObject._who = who
        theObject._tt_cmdH = tt_cmdH
        theObject._nPushPull = n_push_pull
        return theObject

    @staticmethod
    def loadCubeFromIFFMeasureToCreateAn(tt):
        dove = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, tt)
        fits_file_name = os.path.join(dove, 'CubeMeasure.fits')
        hduList = pyfits.open(fits_file_name)
        theObject = AnalyzerIFF()
        cube = np.ma.masked_array(hduList[0].data,
                                  hduList[1].data.astype(bool))
        theObject._cube = cube
        theObject._who = 'Un segmento'
        return theObject

    def setAnalysisMask(self, analysis_mask):
        self._analysisMask = analysis_mask

    def setAnalysisMaskFromMasterMask(self):
        self._analysisMask = self.getMasterMask()

    def setDetectorMask(self, mask_from_ima):
        self._analysisMask = None
        self._rec = None
        self._analysisMask = mask_from_ima

    def getAnalysisMask(self):
        return self._analysisMask

    def _getMaskedInfluenceFunction(self, idx_influence_function):
        return np.ma.array(self.getCube()[:,:,idx_influence_function],
                           mask=self.getAnalysisMask())

    def _createInteractionMatrix(self):
        if self._analysisMask is None:
            self.setAnalysisMaskFromMasterMask()
        n_acts_in_cube = self.getCube().shape[2]
        n_interferometer_pixels_in_mask = \
                    self._getMaskedInfluenceFunction(0).compressed().shape[0]
        self._intMat = np.zeros((n_interferometer_pixels_in_mask,
                                 n_acts_in_cube))
        for i in range(n_acts_in_cube):
            self._intMat[:, i] = \
                            self._getMaskedInfluenceFunction(i).compressed()

    def _createSurfaceReconstructor(self, rCond=1e-15):
        self._rec = self._createRecWithPseudoInverse(rCond)

    def _createRecWithPseudoInverse(self, rCond):
        return np.linalg.pinv(self.getInteractionMatrix(), rcond=rCond)

    def getInteractionMatrix(self):
        if self._intMat is None:
            self._createInteractionMatrix()
        return self._intMat

    def getReconstructor(self):
        if self._rec is None:
            self._createSurfaceReconstructor()
        return self._rec
