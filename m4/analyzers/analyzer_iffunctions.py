'''
Authors
  - C. Selmi: written in 2019
'''

import os
import logging
import h5py
from astropy.io import fits as pyfits
import numpy as np
from m4.ground import read_data
from m4.ground.read_data import InterferometerConverter
from m4.utils.influence_functions_maker import IFFunctionsMaker
from m4.utils.roi import ROI
from m4.utils.img_redux import TipTiltDetrend
from m4.configuration import config_folder_names as fold_name


class AnalyzerIFF():
    '''
    This class analyzes the measurements made through the IFF class by
    generating the cube of measurements and calculating interaction matrix
    and reconstructor.

    HOW TO USE IT::

        from m4.analyzers.analyzer_iffunctions import AnalyzerIFF
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

        self._cmdAmplitude = None
        self._cmdMatrix = None
        self._actsVector = None
        self._nPushPull = None
        self._template = None
        self._h5Folder = None
        self._tt_cmdH = None

    def _ttData(self):
        """ Allow to obtain the tracking number from the path """
        split = os.path.split(self._h5Folder)
        self._tt = split[1]
        return self._tt

    @staticmethod
    def loadInfoFromIFFsTtFolder(tt):
        """ Creates the object using information about path measurements

        Parameters
        ----------
                tt: string
                        measurement tracking number

        Returns
        -------
                theObject: object
                        analyzerIFF class object
        """
        #theObject = AnalyzerIFF()
        theObject = IFFunctionsMaker.loadInfo(tt)
        theObject._cmdAmplitude, theObject._actsVector, theObject._cmdMatrix = \
            read_data.readTypeFromFitsName(theObject._amplitudeTag,
                                           theObject._modesVectorTag,
                                           theObject._cmdMatrixTag) 
        return theObject



    def getCube(self):
        '''
        Returns
        -------
                cube: masked array [pixels, pixels, number of images]
                    cube from analysis
        '''
        return self._cube

    def getMasterMask(self):
        '''
        Returns
        -------
                master_mask: [pixels, pixels]
                            product of the masks of the cube
        '''
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
                where.append(a[0][0])
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
                indvect = where[p]+ indexing_input.shape[0] * k
                vect[indvect] = amplitude[i]
        return vect


    def createCubeFromImageFolder(self, data_file_path=None, tiptilt_detrend=None,
                                  phase_ambiguity=None):
        '''
        Parameters
        ----------
                data_file_path: string
                                measurement data file path

        Other Parameters
        ----------
                ttDetrend: optional
                            in the creation of the cube the images are reduced
                            removing tip tilt on the central segment

                phaseSolve: optional

        Returns
        -------
                cube = masked array [pixels, pixels, number of images]
                        cube from analysis
        '''
        if data_file_path is None:
            data_file_path = self._h5Folder #viene da loadInfoFromIFFsTtFolder
        else:
            data_file_path = data_file_path
        cube_all_act = None
        where = self._indexReorganization()
        ampl_reorg = self._amplitudeReorganization(self._actsVector,
                                                   self._indexingList,
                                                   self._cmdAmplitude,
                                                   self._nPushPull)
        for i in range(self._actsVector.shape[0]):
            print(i)
            for k in range(self._nPushPull):
                p = self._nPushPull * i + k
                n = where[p]
                mis_amp = k* self._indexingList.shape[1] + n
                mis = k * self._indexingList.shape[1] * self._template.shape[0] \
                        + n * self._template.shape[0]

                name = 'img_%04d.h5' %mis
                file_name = os.path.join(data_file_path, name)
                image0 = self._ic.from4D(file_name)

#                 image_sum = np.zeros((image_for_dim.shape[0],
#                                       image_for_dim.shape[1]))
                image_list = [image0]
                for l in range(1, self._template.shape[0]):
                    name = 'img_%04d.h5' %(mis+l)
                    file_name = os.path.join(data_file_path, name)
                    ima = self._ic.from4D(file_name)
                    image_list.append(ima)

                image = np.zeros((image0.shape[0], image0.shape[1]))
                for p in range(1, len(image_list)):
                    #opd2add   = opdtemp[*,*,p]*form[p]+opdtemp[*,*,p-1]*form[p-1]
                    #opd      += opd2add
                    opd2add = image_list[p] * self._template[p] + image_list[p-1] * self._template[p-1]
                    master_mask2add = np.ma.mask_or(image_list[p].mask, image_list[p-1].mask)
                    if p==1:
                        master_mask = master_mask2add
                    else:
                        master_mask = np.ma.mask_or(master_mask, master_mask2add)
                    image += opd2add
                image = np.ma.masked_array(image, mask=master_mask)
                    #image = image + ima * self._template[l] #sbagliato
                img_if = image / (2 * ampl_reorg[mis_amp] * (self._template.shape[0] - 1))
                if tiptilt_detrend is None:
                    img_if = img_if
                else:
                    r = ROI()
                    roi = r.roiGenerator(img_if)
                    tt = TipTiltDetrend()
                    img_if = tt.tipTiltDetrend(img_if, roi, 3)

                if_push_pull_kth = img_if

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


    def _logCubeCreation(self, tiptilt_detrend=None, phase_ambiguity=None):
        """ Use logging function to keep a record of
        how the cube was created """
        who = 'Chissa!'
        if (tiptilt_detrend is None and phase_ambiguity is None):
            self._logger.info('Creation of the IFF cube for %s. Location: %s',
                              who, self._tt)
            self._logger.debug('(Tip and tilt= ignored, \
                                    Phase ambiguity= ignored)')
        elif tiptilt_detrend is None:
            self._logger.info('Creation of the IFF cube for %s. Location: %s',
                              who, self._tt)
            self._logger.debug('(Tip and tilt= ignored, \
                                    Phase ambiguity= resolved)')
        elif phase_ambiguity is None:
            self._logger.info('Creation of the IFF cube for %s. Location: %s',
                              who, self._tt)
            self._logger.debug('(Tip and tilt= removed, \
                                     Phase ambiguity= ignored)')
        else:
            self._logger.info('Creation of the IFF cube for %s. Location: %s',
                              who, self._tt)
            self._logger.debug('(Tip and tilt= removed, \
                                    Phase ambiguity= resolved)')


    def saveCube(self, cube_name, fits_or_h5=0):
        """
        Parameters
        ----------
                cube_name: string
                            name to save the cube
                            example 'Cube.fits'
        """
        tt = self._ttData()
        dove = os.path.join(fold_name.IFFUNCTIONS_ROOT_FOLDER, tt)
        file_name = os.path.join(dove, cube_name)
        if fits_or_h5 == 0:
            header = pyfits.Header()
            header['NPUSHPUL'] = self._nPushPull
            pyfits.writeto(file_name, self._cube.data, header)
            pyfits.append(file_name, self._cube.mask.astype(int))
            pyfits.append(file_name, self._cmdAmplitude)
            pyfits.append(file_name, self._actsVector)
        else:
            hf = h5py.File(file_name, 'w')
            hf.create_dataset('dataset_1', data=self._cube.data)
            hf.create_dataset('dataset_2', data=self._cube.mask.astype(int))
            hf.create_dataset('dataset_3', data=self._cmdAmplitude)
            hf.create_dataset('dataset_4', data=self._actsVector)
            hf.attrs['NPUSHPUL'] = self._nPushPull
            hf.close()


    @staticmethod
    def loadAnalyzer(file_name, fits_or_h5=0):
        """ Creates the object using information contained in Cube

        Parameters
        ----------
                fits_file_name: string
                                cube file name path

        Returns
        -------
                theObject: object
                            analyzerIFF class object
        """
        theObject = AnalyzerIFF()
        if fits_or_h5 == 0:
            header = pyfits.getheader(file_name)
            hduList = pyfits.open(file_name)
            cube = np.ma.masked_array(hduList[0].data,
                                      hduList[1].data.astype(bool))
            acts_vector = hduList[3].data
            cmd_amplitude = hduList[2].data
            try:
                n_push_pull = header['NPUSHPUL']
            except KeyError:
                n_push_pull = 1

            theObject._actsVector = acts_vector
            theObject._cmdAmplitude = cmd_amplitude
            theObject._nPushPull = n_push_pull
            theObject._cube = cube
        else:
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
            hf.close()
        return theObject

    def setAnalysisMask(self, analysis_mask):
        ''' Set the analysis mask chosen

        Parameters
        ----------
                analysis_mask: numpy array [pixels, pixels]
        '''
        self._analysisMask = analysis_mask

    def setAnalysisMaskFromMasterMask(self):
        ''' Set the analysis mask using the master mask of analysis cube
        '''
        self._analysisMask = self.getMasterMask()

    def setDetectorMask(self, mask_from_ima):
        ''' Set the detector mask chosen

        Parameters
        ----------
        detector_mask: numpy array [pixels, pixels]
        '''
        self._analysisMask = None
        self._rec = None
        self._analysisMask = mask_from_ima

    def getAnalysisMask(self):
        '''
        Returns
        -------
            analysis_mask: numpy array [pixels, pixels]
        '''
        return self._analysisMask

    def _getMaskedInfluenceFunction(self, idx_influence_function):
        return np.ma.array(self.getCube()[:, :, idx_influence_function],
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
        '''
        Returns
        -------
                intMat: numpy array
                        interaction matrix from cube
        '''
        if self._intMat is None:
            self._createInteractionMatrix()
        return self._intMat

    def getReconstructor(self):
        '''
        Returns
        -------
                rec = numpy array
                    reconstructor calculated as pseudo inverse of the interaction matrix
        '''
        if self._rec is None:
            self._createSurfaceReconstructor()
        return self._rec
