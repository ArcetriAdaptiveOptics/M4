'''
Authors
  - C. Selmi:  written in 2019
'''
import os
import logging
from astropy.io import fits as pyfits
import numpy as np
from m4.configuration.config import fold_name
from m4.configuration.ott_parameters import OttParameters
from m4.ground import tracking_number_folder
from m4.ground import zernike


class OpticalCalibration():
    """
    Class for the optical calibration

    HOW TO USE IT::

        from m4.utils.optical_calibration import OpticalCalibration
        cal = OpticalCalibration()
    """

    def __init__(self, ott, interf):
        """The constructor """
        self._logger = logging.getLogger('OPT_CALIB:')
        self._interf = interf
        self._ott = ott
        self._cube = None
        self._rec = None
        self._intMat = None
        self._intMatNorm = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return fold_name.CALIBRATION_ROOT_FOLDER


    def measureCalibrationMatrix(self, who, command_amp_vector, n_push_pull,
                                 n_frames, old_or_new):
        '''
        Parameters
        ----------
                ott = object
                    tower
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
                n_frames: int
                        number of frame for 4D measurement
                old_or_new: int
                        0 for new (mixed), 1 for old (not mixed)

        Returns
        -------
        tt : string
            tracking number
        '''
        #self._ott = ott
        self._nPushPull = n_push_pull
        self._commandAmpVector = command_amp_vector
        self._who = who
        self._old_or_new = old_or_new

        save = tracking_number_folder.TtFolder(self._storageFolder())
        dove, self._tt = save._createFolderToStoreMeasurements()
        self._logger.info('Measure of calibration. Location: %s', self._tt)
        self._saveCommandAmplitudeAsFits(dove)

        self._commandMatrix, self._commandList = self._createCommandMatrix(who,
                                                                           command_amp_vector,
                                                                           old_or_new)
        self._saveCommandMatrixAsFits(dove)
        self._saveCMat(dove)
        self._measureAndStore(who, self._commandList, n_push_pull, dove, n_frames)
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
        a = OpticalCalibration.loadCommandMatrixFromFits(tt)
        #cubo non normalizzato
        a.createCube(tt, 0)
        cube = a.getCube()

        ima = cube[:,:,0]
        from m4.utils.roi import ROI
        r = ROI()
        roi = r.roiGenerator(ima)
        mask = roi[mask_index]
        dove = os.path.join(self._storageFolder(), tt)
        self._saveMask(dove, mask)

        self._intMat = a.getInteractionMatrix(mask)
        fits_file_name = os.path.join(dove, 'InteractionMatrix.fits')
        pyfits.writeto(fits_file_name, self._intMat, overwrite=True)

        #cubo normalizzato
        if a._old_or_new == 1:
            a.createCube(tt)
            cube_norm = a.getCube()
            self._intMatNorm = a.getInteractionMatrix(mask)
            fits_file_name = os.path.join(dove, 'InteractionMatrixNorm.fits')
            pyfits.writeto(fits_file_name, self._intMatNorm, overwrite=True)
        return self._intMat, self._intMatNorm

    def _measureAndStore(self, who, command_list, n_push_pull, dove, n_frames):
        if who == 0:
            vec_push_pull = np.array((1, -1))
            #mis = (len(command_list)-2) * n_push_pull * vec_push_pull.size
            par0 = self._ott.parabola.getPosition()
            rm0 = self._ott.referenceMirror.getPosition()
            for k in range(n_push_pull):
                for i in range(len(command_list)-2):
                    j = (len(command_list)-2)*k *2
                    mis = np.array([j , j +1])
                    if i==0:
                        pcmd = np.array(command_list[i])
                        for v in range(vec_push_pull.size):
                            par1 = pcmd * vec_push_pull[v]
                            print(par1)
                            self._ott.parabola.setPosition(par0+par1)
                            masked_ima = self._interf.acquire_phasemap(n_frames)
                            name = 'Frame_%04d.fits' %(2*i+mis[v])
                            print(name)
                            self._interf.save_phasemap(dove, name, masked_ima)
                            self._ott.parabola.setPosition(par0)
                    elif i==1 or i==2:
                        if i==1:
                            l=i
                        else:
                            l=i+1
                        pcmd = np.array(command_list[l])
                        rcmd = np.array(command_list[l+1])
                        for v in range(vec_push_pull.size):
                            par1 = pcmd * vec_push_pull[v]
                            rm1 = rcmd * vec_push_pull[v]
                            self._ott.parabola.setPosition(par0+par1)
                            if np.count_nonzero(rm1) !=0:
                                self._ott.referenceMirror.setPosition(rm0+rm1)
                            print(par1, rm1)
                            masked_ima = self._interf.acquire_phasemap(n_frames)
                            name = 'Frame_%04d.fits' %(2*i+mis[v])
                            print(name)
                            self._interf.save_phasemap(dove, name, masked_ima)
                            self._ott.parabola.setPosition(par0)
                            self._ott.referenceMirror.setPosition(rm0)
                    else:
                        rcmd = np.array(command_list[i+2])
                        for v in range(vec_push_pull.size):
                            rm1 = rcmd * vec_push_pull[v]
                            self._ott.referenceMirror.setPosition(rm0+rm1)
                            print(rm1)
                            masked_ima = self._interf.acquire_phasemap(n_frames)
                            name = 'Frame_%04d.fits' %(2*i+mis[v])
                            print(name)
                            self._interf.save_phasemap(dove, name, masked_ima)
                            self._ott.referenceMirror.setPosition(rm0)
        elif who == 1:
            pass
        elif who == 2:
            pass
        elif who == 3:
            #m4 calib
            pass


    def _createCommandMatrix(self, who, command_amp_vector, old_or_new):
        '''
            args:
                who=
                    0 per mixing
                    1 per parable
                    2 per reference mirror
                    3 per deformable mirror
                command_amp_vector = vector of command amplitude
                old_or_new: int
                        0 for new (mixed), 1 for old (not mixed)

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
        commandMatrix, commandList = self._createCmatAndCmdList(command_amp_vector,
                                                                self._dofIndex, old_or_new)

        return commandMatrix, commandList

    def _createCmatAndCmdList(self, command_amp_vector, dofIndex_vector, old_or_new):
        '''
        Parameters
        ---------
        old_or_new = 0 for new (mixed), 1 for old (not mixed)
        '''
        if old_or_new == 0:
            #crea matrice 5 x 5
            command_matrix = np.zeros((command_amp_vector.size, command_amp_vector.size))
            for i in range(command_amp_vector.shape[0]):
                j = i
                if i==1 or i==2:
                    command_matrix[i,j] = command_amp_vector[i]
                    command_matrix[i,j+2] = OttParameters.par_rm_coef_for_coma_measuremets*command_amp_vector[i]
                else:
                    command_matrix[i,j] = command_amp_vector[i]
        elif old_or_new == 1:
            command_matrix = np.zeros((command_amp_vector.size, command_amp_vector.size))
            for i in range(command_amp_vector.shape[0]):
                command_matrix[i,i] = command_amp_vector[i]

        #crea i comandi
        command_list = []
        for i in range(command_matrix.shape[0]):
            if i==1 or i==2:
                cmd = np.zeros(6)
                cmd[dofIndex_vector[i]] = command_matrix[i,i]
                command_list.append(cmd)
                cmd1 = np.zeros(6)
                cmd1[dofIndex_vector[i+2]] = command_matrix[i,i+2]
                command_list.append(cmd1)
            else:
                cmd_amp = command_matrix[i,i]
                cmd = np.zeros(6)
                cmd[dofIndex_vector[i]] = command_matrix[i,i]
                command_list.append(cmd)
        return command_matrix, command_list

    def _saveCommandAmplitudeAsFits(self, dove):
        fits_file_name = os.path.join(dove, 'CommandAmplitude.fits')
        pyfits.writeto(fits_file_name, self._commandAmpVector)

    def _saveCMat(self, dove):
        fits_file_name = os.path.join(dove, 'CMat.fits')
        pyfits.writeto(fits_file_name, self._commandMatrix.T)

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
        header['OLDNEW'] = self._old_or_new
        pyfits.writeto(fits_file_name, self._commandAmpVector, header)
        pyfits.append(fits_file_name, self._commandMatrix.T, header)

    @staticmethod
    def loadCommandMatrixFromFits(tt):
        """ Creates the object using information contained in command matrix fits file

        Parameters
        ----------
                tt: string
                    tracking number

        Returns
        -------
                theObject: ibjecct
                         opt_calibration class object
        """
        ott = None
        interf = None
        theObject = OpticalCalibration(ott, interf)
        theObject._tt = tt
        dove = os.path.join(theObject._storageFolder(), tt)
        file = os.path.join(dove, 'CommandMatrixInfo.fits')
        header = pyfits.getheader(file)
        hduList = pyfits.open(file)
        theObject._who = header['WHO']
        theObject._nPushPull = header['NPUSHPUL']
        theObject._old_or_new = header['OLDNEW']
        theObject._commandAmpVector = hduList[0].data
        theObject._commandMatrix = hduList[1].data
        return theObject

    def createCube(self, tt, norm=1):
        """
        Parameters
        ----------
            tt = tracking number
        """
        if norm != 1:
            self._commandAmpVector = np.ones(self._commandAmpVector.size)
        self._logger.info('Creation of the cube relative to %s', tt)
        self._cube = None
        fold = os.path.join(OpticalCalibration._storageFolder(), tt)
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
        if self._cube is None:
            self.createCube(self._tt)
            self._cube = self.getCube()
        return self._cube

    def _createInteractionMatrixAndReconstructor(self, mask):
        coefList = []
        for i in range(self._cube.shape[2]):
            ima = np.ma.masked_array(self._cube[:,:,i], mask=mask)
            coef, mat = zernike.zernikeFit(ima, np.arange(10)+1)
            #z= np.array([2,3,4,7,8])
            z = np.array([1, 2, 3, 6, 7])
            final_coef = np.zeros(z.shape[0])
            final_coef = coef[z]
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

    def _saveMask(self, dove, mask):
        """
        args:
            dove = path that indicates where to save the mask file
        """
        fits_file_name = os.path.join(dove, 'Mask.fits')
        pyfits.writeto(fits_file_name, mask.astype(int), overwrite=True)
