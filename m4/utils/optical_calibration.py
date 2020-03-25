'''
@author: cselmi
'''
import os
import logging
from astropy.io import fits as pyfits
import numpy as np
from m4.ground.configuration import Configuration
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from m4.ground import tracking_number_folder


class opt_calibration():
    """
    Class for the optical calibration

    HOW TO USE IT:
    from m4.utils.optical_calibration import opt_calibration
    cal = opt_calibration()
    """

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('OPT_CALIB:')
        self._zOnM4 = ZernikeOnM4()
        self._rec = None
        self._intMat = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return os.path.join(Configuration.OPD_DATA_FOLDER,
                            "Calibration")


    def measureCalibrationMatrix(self, who, command_amp_vector, n_push_pull):
        '''
            arg:
                who = number indicating the optical element
                    on which to perform the calibration
                    0 per mixing
                    1 per parable
                    2 per reference mirror
                    3 per deformable mirror
                command_amp_vector = vector containing the amplitude of the
                                    commands to give degrees of freedom to
                                    calibrate
                n_push_pull = number of push pull

            Returns:
                    tt = tracking number
        '''
        self._nPushPull = n_push_pull
        self._commandAmpVector = command_amp_vector

        store_in_folder = self._storageFolder()
        save = tracking_number_folder.TtFolder(store_in_folder)
        dove, self._tt = save._createFolderToStoreMeasurements()
        self._logger.info('Measure of calibration. Location: %s', self._tt)

        self._commandMatrix = self._createCommandMatrix(who,
                                                        self._commandAmpVector,
                                                        self._nPushPull)
        self._saveCommandMatrixAsFits(dove)

        self._applyCommandMatrix(self._commandMatrix)
        return self._tt

    def analyzerCalibrationMeasurement(self, tt, mask_index):
        '''
        args:
            tt = tracking number of the measures to be analysed
            mask_index = int the reference mirror mask index

        returns:
                self_intMat = interation matrix
                self._rec = reconstructor
        '''
        a = opt_calibration.loadCommandMatrixFromFits(tt)
        a.createCube(tt)
        cube = a.getCube()
        ima = cube[:,:,0]
        from m4.utils.roi import ROI
        r = ROI()
        roi = r.roiGenerator(ima)
        mask = roi[mask_index]
        dove = os.path.join(self._storageFolder(), tt)
        self._saveMask(dove, mask)
        self._intMat = a.getInteractionMatrix(mask)
        self._rec = a.getReconstructor(mask)
        self._saveIntMatAndRec(dove)
        return self._intMat, self._rec


    def _applyCommandMatrix(self, command_matrix):
        #deve applicare la matrice e salvare gli interferogrammi
        pass


    def _createCommandMatrix(self, who, command_amp_vector, n_push_pull):
        '''
            args:
                who=
                    0 per mixing
                    1 per parable
                    2 per reference mirror
                    3 per deformable mirror
                command_amp_vector = vector of command amplitude
                n_push_pull = number of push pull

            returns:
                    self._commandMatrix = command matrix for the dof
        '''
        if who == 0:
            self._who = 'PAR + RM'
            self._rows = Configuration.PARABOLA_DOF + Configuration.RM_DOF
        elif who == 1:
            self._who = 'PAR'
            self._rows = Configuration.PARABOLA_DOF
        elif who == 2:
            self._who = 'RM'
            self._rows = Configuration.RM_DOF
        elif who == 3:
            self._who = 'M4'
            self._rows = Configuration.M4_DOF
        else:
            raise OSError('Who= %s doesnt exists' % who)

        self._logger.info('Creation of the command matrix for %s', self._who)
        self._commandMatrix = self._createCommandHistoryMatrix(self._rows,
                                                               command_amp_vector,
                                                               n_push_pull)
        return self._commandMatrix


    def _createCommandHistoryMatrix(self, rows, command_amp_vector,
                                    n_push_pull):
        '''
            create the command matrix using as lines
            the degrees of freedom chosen in rows

        args:
            rows = int
            command_amp_vector = vector of command amplitude
            n_push_pull = number of push pull

            returns:
                    command_matrix = command matrix for the dof
        '''
        vec_push_pull = np.array((1, -1))
        rows = rows
        columns = command_amp_vector.shape[0]*n_push_pull*vec_push_pull.shape[0]
        command_matrix = np.zeros((rows, columns))

        for k in range(n_push_pull):
            for i in range(command_amp_vector.shape[0]):
                j = (command_amp_vector.shape[0]*vec_push_pull.shape[0])*k +2*i
                command_matrix[i,j] = command_amp_vector[i]*vec_push_pull[0]
                command_matrix[i,j+1] = command_amp_vector[i]*vec_push_pull[1]
        return command_matrix


    def _saveCommandMatrixAsFits(self, dove):
        """
        Save fits file for the command matrix and the data relating to
        its creation

        args:
            dove = path that indicates where to save the command matrix file
        """
        fits_file_name = os.path.join(dove, 'CommandMatrix.fits')
        header = pyfits.Header()
        header['NPUSHPUL'] = self._nPushPull
        header['WHO'] = self._who
        pyfits.writeto(fits_file_name, self._commandAmpVector, header)
        pyfits.append(fits_file_name, self._commandMatrix, header)

    @staticmethod
    def loadCommandMatrixFromFits(tt):
        """ Creates the object using information contained in command matrix
            fits file

            Args:
                tt = tracking number

            Return:
                theObject = opt_calibration class object
        """
        theObject = opt_calibration()
        theObject._tt = tt
        dove = os.path.join(theObject._storageFolder(), tt)
        file = os.path.join(dove, 'CommandMatrix.fits')
        header = pyfits.getheader(file)
        hduList = pyfits.open(file)
        theObject._who = header['WHO']
        theObject._nPushPull = header['NPUSHPUL']
        theObject._commandAmpVector = hduList[0].data
        theObject._commandMatrix = hduList[1].data
        return theObject

    def _testCalibration_createCubeMeasurefromFileFitsMeasure(self):
        """
        Test function for the cube measure creation
        """
        cube_measure = None
        fold_my_pc = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/MixingIntMat/20190930_162714'
        #fold_m4_pc = os.path.join(Configuration.OPD_DATA_FOLDER, 'TestData', 'MixingIntMat')
        for i in range(5):
            name = 'Frame_%04d.fits' %i
            file = os.path.join(fold_my_pc, name)
            hduList = pyfits.open(file)
            aa = hduList[0].data
            mode = np.ma.masked_array(aa[0,:,:], mask=np.invert(aa[1,:,:].astype(bool)))
            if cube_measure is None:
                cube_measure = mode
            else:
                cube_measure = np.ma.dstack((cube_measure, mode))
        return cube_measure

    def createCube(self, tt):
        """
        args:
            tt = tracking number
        """
        self._logger.info('Creation of the cube relative to %s', tt)
        cube_from_measure = self._testCalibration_createCubeMeasurefromFileFitsMeasure()
        for i in range(cube_from_measure.shape[2]):
            cube_from_measure[:,:,i] = cube_from_measure[:,:,i] / self._commandAmpVector[i]
        self._cube = cube_from_measure
        return

    def getCube(self):
        return self._cube

    def _createInteractionMatrixAndReconstructor(self, mask):
        coefList = []
        for i in range(self._cube.shape[2]):
            ima = np.ma.masked_array(self._cube[:,:,i], mask=mask)
            coef, mat = self._zOnM4.zernikeFit(ima, np.arange(2, 11))
            #z= np.array([2,3,4,7,8])
            z = np.array([0, 1, 2, 5, 6])
            final_coef = np.zeros(z.shape[0])
            aa = np.arange(final_coef.shape[0])
            zipped = zip(aa, z)
            for i, j in zipped:
                final_coef[i] = coef[j]
            self._mat = mat
            coefList.append(final_coef)

        self._intMat = np.zeros((coefList[0].shape[0], self._cube.shape[2]))
        for j in range(self._cube.shape[2]):
            self._intMat.T[j] = coefList[j]

        self._rec = np.linalg.pinv(self._intMat)


    def getInteractionMatrix(self, mask):
        if self._intMat is None:
            self._createInteractionMatrixAndReconstructor(mask)
        return self._intMat

    def getReconstructor(self, mask):
        if self._rec is None:
            self._createInteractionMatrixAndReconstructor(mask)
        return self._rec

    def _saveIntMatAndRec(self, dove):
        """
        args:
            dove = path that indicates where to save the files
        """
        fits_file_name = os.path.join(dove, 'InteractionMatrix.fits')
        pyfits.writeto(fits_file_name, self._intMat)
        fits_file_name = os.path.join(dove, 'Reconstructor.fits')
        pyfits.writeto(fits_file_name, self._rec)

    def _saveMask(self, dove, mask):
        """
        args:
            dove = path that indicates where to save the mask file
        """
        fits_file_name = os.path.join(dove, 'Mask.fits')
        pyfits.writeto(fits_file_name, mask.astype(int))
