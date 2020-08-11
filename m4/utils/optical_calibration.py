'''
@author: cselmi
'''
import os
import logging
from astropy.io import fits as pyfits
import numpy as np
from m4.configuration.config import fold_name
from m4.configuration.ott_parameters import OttParameters
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from m4.ground import tracking_number_folder
from m4.utils.interface_4D import comm4d


class opt_calibration():
    """
    Class for the optical calibration

    HOW TO USE IT::

        from m4.utils.optical_calibration import opt_calibration
        cal = opt_calibration()
    """

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('OPT_CALIB:')
        self._zOnM4 = ZernikeOnM4()
        self._c4d = comm4d()
        self._rec = None
        self._intMat = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return fold_name.CALIBRATION_ROOT_FOLDER


    def measureCalibrationMatrix(self, ott, who, command_amp_vector, n_push_pull):
        '''
        Parameters
        ----------
                who: int
                     number indicating the optical element
                    on which to perform the calibration
                    0 for mixing
                    1 for parable
                    2 for reference mirror
                    3 for deformable mirror
                command_amp_vector: numpy array
                                    vector containing the amplitude of the
                                    commands to give degrees of freedom to
                                    calibrate
                n_push_pull: int
                            number of push pull

        Returns
        -------
        tt : string
            tracking number
        '''
        self._ott = ott
        self._nPushPull = n_push_pull
        self._commandAmpVector = command_amp_vector
        self._who = who

        store_in_folder = self._storageFolder()
        save = tracking_number_folder.TtFolder(store_in_folder)
        dove, self._tt = save._createFolderToStoreMeasurements()
        self._logger.info('Measure of calibration. Location: %s', self._tt)
        self._saveCommandAmplitudeAsFits(dove)

        self._commandMatrix = self._createCommandMatrix(who,
                                                        self._commandAmpVector,
                                                        self._nPushPull)
        self._saveCommandMatrixAsFits(dove)

        self._measureAndStoreCommandMatrix(who, self._commandMatrix, dove)
        return self._tt

    def analyzerCalibrationMeasurement(self, tt, mask_index):
        '''
        Parameters
        ----------
            tt: string
                tracking number of the measures to be analysed
            mask_index: int
                    reference mirror mask index

        Returns
        -------
                self_intMat: numpy array
                         interation matrix
                self._rec: numpy array
                         reconstructor
        '''
        #mask_index per il simulatore: 2 (se le ruoto: 3)
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


    def _measureAndStoreCommandMatrix(self, who, command_matrix, dove):
        #deve applicare la matrice e salvare gli interferogrammi
        command_list = []
        for i in range(command_matrix.shape[1]):
            cmd = command_matrix[:,i]
            command_list.append(cmd)
        if who == 0:
#             self._ott.slide(0.75)
#             self._ott.angle(90.)
#             self._ott.rslide(0.6)
            for l in range(self._nPushPull):
                for m in range(np.int(len(command_list)/self._nPushPull)):
                    k = 2 * l * self._dofIndex.size + m
                    if k <= 2 * l * self._dofIndex.size + self._dofIndex.size:
                        self._ott.parab(command_list[k])
                    elif 2 * l * self._dofIndex.size + self._dofIndex.size < k < 2 * (l+1) * self._dofIndex.size:
                        self._ott.parab(np.zeros(6))
                        self._ott.refflat(command_list[k])
                    masked_ima = self._c4d.acq4d(self._ott, 1, show=1)
                    #masked_ima = np.ma.masked_array(p, mask=np.invert(m.astype(bool)))
                    name = 'Frame_%04d.fits' %k
                    self._saveSimulatedInterf(dove, name, masked_ima)
                self._ott.refflat(np.zeros(6))

        elif who == 1:
            pass
        elif who == 2:
            pass
        elif who == 3:
#             self._ott.slide(0.75)
#             self._ott.angle(90.)
#             self._ott.rslide(0.6)
            for k in range(len(command_list)):
                self._ott.m4(-command_list[i])
                masked_ima = self._c4d.acq4d(self._ott, 1, show=1)
                name = 'Frame_%04d.fits' %k
                self._saveSimulatedInterf(dove, name, masked_ima)
            self._ott.m4(np.zeros(6))

    def _saveSimulatedInterf(self, dove, file_name, image):
        fits_file_name = os.path.join(dove, file_name)
        pyfits.writeto(fits_file_name, image.data)
        pyfits.append(fits_file_name, image.mask.astype(int))


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
            self._dofIndex = np.append(OttParameters.PARABOLA_DOF, OttParameters.RM_DOF)
        elif who == 1:
            self._who = 'PAR'
            self._dofIndex = OttParameters.PARABOLA_DOF
        elif who == 2:
            self._who = 'RM'
            self._dofIndex = OttParameters.RM_DOF
        elif who == 3:
            self._who = 'M4'
            self._dofIndex = OttParameters.M4_DOF
        else:
            raise OSError('Who= %s doesnt exists' % who)

        self._logger.info('Creation of the command matrix for %s', self._who)
        self._commandMatrix = self._createCommandHistoryMatrix(self._dofIndex,
                                                               command_amp_vector,
                                                               n_push_pull)
        return self._commandMatrix


    def _createCommandHistoryMatrix(self, dofIndex_vector, command_amp_vector, n_push_pull):
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
        rows = 6
        columns = command_amp_vector.shape[0]*n_push_pull*vec_push_pull.shape[0]
        command_matrix = np.zeros((rows, columns))

        for k in range(n_push_pull):
            for i in range(command_amp_vector.shape[0]):
                j = (command_amp_vector.shape[0]*vec_push_pull.shape[0])*k +2*i
                command_matrix[dofIndex_vector[i],j] = command_amp_vector[i]*vec_push_pull[0]
                command_matrix[dofIndex_vector[i],j+1] = command_amp_vector[i]*vec_push_pull[1]
        return command_matrix

    def _saveCommandAmplitudeAsFits(self, dove):
        fits_file_name = os.path.join(dove, 'CommandAmplitude.fits')
        pyfits.writeto(fits_file_name, self._commandAmpVector)

    def _saveCommandMatrixAsFits(self, dove):
        """
        Save fits file for the command matrix and the data relating to
        its creation

        args:
            dove = path that indicates where to save the command matrix file
        """
        fits_file_name = os.path.join(dove, 'CommandMatrixInfo.fits')
        header = pyfits.Header()
        header['NPUSHPUL'] = self._nPushPull
        header['WHO'] = self._who
        pyfits.writeto(fits_file_name, self._commandAmpVector, header)
        pyfits.append(fits_file_name, self._commandMatrix, header)

    @staticmethod
    def loadCommandMatrixFromFits(tt):
        """ Creates the object using information contained in command matrix
            fits file

        Parameters
        ----------
                tt: string
                    tracking number

        Returns
        -------
                theObject: ibjecct
                         opt_calibration class object
        """
        theObject = opt_calibration()
        theObject._tt = tt
        dove = os.path.join(theObject._storageFolder(), tt)
        file = os.path.join(dove, 'CommandMatrixInfo.fits')
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
        amp = np.array([5.0e-06, 5.0e-06, 2.5e-05, 5.0e-06, 5.0e-06, 5.0e-06])
        for i in range(5):
            name = 'Frame_%04d.fits' %i
            file = os.path.join(fold_my_pc, name)
            hduList = pyfits.open(file)
            aa = hduList[0].data
            mode = np.ma.masked_array(aa[0,:,:], mask=np.invert(aa[1,:,:].astype(bool)))
            mode_norm = mode/amp[i]
            if cube_measure is None:
                cube_measure = mode_norm
            else:
                cube_measure = np.ma.dstack((cube_measure, mode_norm))
        return cube_measure

    def createCube(self, tt):
        """
        Parameters
        ----------
            tt = tracking number
        """
        self._logger.info('Creation of the cube relative to %s', tt)
#         cube_from_measure = self._testCalibration_createCubeMeasurefromFileFitsMeasure()
#         for i in range(cube_from_measure.shape[2]):
#             cube_from_measure[:,:,i] = cube_from_measure[:,:,i] / self._commandAmpVector[i]
#         self._cube = cube_from_measure
        self._cube = None
        fold = os.path.join(opt_calibration._storageFolder(), tt)
        for i in range(self._commandAmpVector.shape[0]):
            for j in range(self._nPushPull):
                k = 2*i + 2*self._commandAmpVector.shape[0]*j
                name_pos = 'Frame_%04d.fits' %k
                name_neg = 'Frame_%04d.fits' %(k+1)
                file = os.path.join(fold, name_pos)
                hduList = pyfits.open(file)
                image_pos = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
                file = os.path.join(fold, name_neg)
                hduList = pyfits.open(file)
                image_neg = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))

                if self._who == 'PAR + RM':
                    image = (image_pos - image_neg) / (2 * self._commandAmpVector[i])
                elif self._who == 'M4':
                    image = (image_pos + image_neg) / (2 * self._commandAmpVector[i])

                if j == 0:
                    all_push_pull_act = image
                else:
                    all_push_pull_act = np.ma.dstack((all_push_pull_act, image))
            if self._nPushPull == 1:
                final_ima = all_push_pull_act
            else:
                final_ima = np.ma.mean(all_push_pull_act, axis=2)

            if self._cube is None:
                self._cube = final_ima
            else:
                self._cube = np.ma.dstack((self._cube, final_ima))
        return

    def getCube(self):
        '''
        Returns
        -------
        cube: numpy masked array
            analyzed measurements
        '''
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
        '''
        Parameters
        ----------
        mask: numpy array
            reference mirror mask

        Returns
        -------
        intMat: numpy array
                interaction matrix
        '''
        if self._intMat is None:
            self._createInteractionMatrixAndReconstructor(mask)
        return self._intMat

    def getReconstructor(self, mask):
        '''
        Parameters
        ----------
        mask: numpy array
            reference mirror mask

        Returns
        -------
        rec: numpy array
            reconstructor
        '''
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
