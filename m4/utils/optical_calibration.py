'''
Authors
  - C. Selmi:  written in 2019
'''
import os
import logging
import time
from astropy.io import fits as pyfits
import numpy as np
from m4.configuration.config import fold_name
from m4.configuration.ott_parameters import OttParameters
from m4.ground import tracking_number_folder
from m4.ground.interface_4D import comm4d
from m4.ground import zernike


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
        self._c4d = comm4d()
        self._cube = None
        self._rec = None
        self._intMat = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return fold_name.CALIBRATION_ROOT_FOLDER


    def measureCalibrationMatrix(self, ott, who, command_amp_vector, n_push_pull, n_frames, old=1):
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

        if old==1:
            self._commandMatrix = self._createCommandMatrix(who,
                                                            self._commandAmpVector,
                                                            self._nPushPull)
            self._saveCommandMatrixAsFits(dove)
            self._measureAndStoreCommandMatrix(who, self._commandMatrix, dove, n_frames)
        else:
            self._commandMatrix, self._commandList = self._createCommandMatrix(who,
                                                            self._commandAmpVector,
                                                            self._nPushPull, old)
            self._saveCommandMatrixAsFits(dove)
            self._measureAndStoreMix(self._commandList, self._nPushPull, dove, n_frames)
        return self._tt

    def analyzerCalibrationMeasurement(self, tt, mask_index, norm=1):
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
        a.createCube(tt, norm)
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

    def _measureAndStoreMix(self, command_list, n_push_pull, dove, n_frames):
        vec_push_pull = np.array((1, -1))
        #mis = (len(command_list)-2) * n_push_pull * vec_push_pull.size
        par0 = self._ott.parab()
        rm0 = self._ott.refflat()
        for k in range(n_push_pull):
            for i in range(len(command_list)-2):
                j = (len(command_list)-2)*k *2
                mis = np.array([j , j +1])
                if i==0:
                    pcmd = np.array(command_list[i])
                    for v in range(vec_push_pull.size):
                        par1 = pcmd * vec_push_pull[v]
                        print(par1)
                        self._ott.parab(par0+par1)
                        masked_ima = self._c4d.acq4d(n_frames, 0, self._ott)
                        name = 'Frame_%04d.fits' %(2*i+mis[v])
                        print(name)
                        self._c4d.save_phasemap(dove, name, masked_ima)
                        self._ott.parab(par0)
                elif i==1 or i==2:
                    pcmd = np.array(command_list[i])
                    rcmd = np.array(command_list[i+1])
                    for v in range(vec_push_pull.size):
                        par1 = pcmd * vec_push_pull[v]
                        rm1 = rcmd * vec_push_pull[v]
                        self._ott.parab(par0+par1)
                        self._ott.refflat(rm0+rm1)
                        print(par1, rm1)
                        masked_ima = self._c4d.acq4d(n_frames, 0, self._ott)
                        name = 'Frame_%04d.fits' %(2*i+mis[v])
                        print(name)
                        self._c4d.save_phasemap(dove, name, masked_ima)
                        self._ott.parab(par0)
                        self._ott.refflat(rm0)
                else:
                    rcmd = np.array(command_list[i+2])
                    for v in range(vec_push_pull.size):
                        rm1 = rcmd * vec_push_pull[v]
                        self._ott.refflat(rm0+rm1)
                        print(rm1)
                        masked_ima = self._c4d.acq4d(n_frames, 0, self._ott)
                        name = 'Frame_%04d.fits' %(2*i+mis[v])
                        print(name)
                        self._c4d.save_phasemap(dove, name, masked_ima)
                        self._ott.refflat(rm0)


    def _measureAndStoreCommandMatrix(self, who, command_matrix, dove, n_frames):
        #deve applicare la matrice e salvare gli interferogrammi
        command_list = []
        for i in range(command_matrix.shape[1]):
            cmd = command_matrix[:,i]
            command_list.append(cmd)
        if who == 0:
            par0 = self._ott.parab()
            rm0 = self._ott.refflat()
            for l in range(self._nPushPull):
                for m in range(np.int(len(command_list)/self._nPushPull)):
                    k = 2 * l * self._dofIndex.size + m
                    if k <= 2 * l * self._dofIndex.size + self._dofIndex.size:
                        self._ott.parab(par0 + command_list[k])
                    elif 2 * l * self._dofIndex.size + self._dofIndex.size < k < 2 * (l+1) * self._dofIndex.size:
                        self._ott.parab(par0)
                        self._ott.refflat(rm0 + command_list[k])
                    #time.sleep(5)
                    #print('5 secondi di attesa')
                    masked_ima = self._c4d.acq4d(n_frames, self._ott)
                    #masked_ima = np.ma.masked_array(p, mask=np.invert(m.astype(bool)))
                    name = 'Frame_%04d.fits' %k
                    self._c4d.save_phasemap(dove, name, masked_ima)
                self._ott.refflat(rm0)

        elif who == 1:
            pass
        elif who == 2:
            pass
        elif who == 3:
            m40 = self._ott.m4()
            for k in range(len(command_list)):
                self._ott.m4(m40-command_list[i])
                masked_ima = self._c4d.acq4d(1, self._ott)
                name = 'Frame_%04d.fits' %k
                self._c4d.save_phasemap(dove, name, masked_ima)
            self._ott.m4(m40)


    def _createCommandMatrix(self, who, command_amp_vector, n_push_pull, old=1):
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
        if old==1:
            self._commandMatrix = self._createCommandHistoryMatrix(self._dofIndex,
                                                                   command_amp_vector,
                                                                   n_push_pull)
            return self._commandMatrix
        else:
            self._commandMatrix, self._commandList = self._createCommandMatrixMix(command_amp_vector, self._dofIndex)
            return self._commandMatrix, self._commandList

    def _createCommandMatrixMix(self, command_amp_vector, dofIndex_vector):
        #crea matrice 5 x 5
        command_matrix = np.zeros((command_amp_vector.size, command_amp_vector.size))
        for i in range(command_amp_vector.shape[0]):
            j = i
            if i==1 or i==2:
                command_matrix[i,j] = command_amp_vector[i]
                command_matrix[i,j+2] = -2.5*command_amp_vector[i]
            else:
                command_matrix[i,j] = command_amp_vector[i]
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
#         par_cmd = [command_list[0], command_list[1], command_list[3]]
#         rm_cmd = [command_list[2], command_list[4], command_list[5], command_list[6]]  
        return command_matrix, command_list

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

    def createCube(self, tt, norm=1):
        """
        Parameters
        ----------
            tt = tracking number
        """
        if norm !=1:
            self._commandAmpVector = np.ones(self._commandAmpVector.size)
        self._logger.info('Creation of the cube relative to %s', tt)
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
        if self._cube is None:
            self._cube = self.createCube(self._tt)
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

    def _saveIntMatAndRec(self, dove):
        """
        args:
            dove = path that indicates where to save the files
        """
        fits_file_name = os.path.join(dove, 'InteractionMatrix.fits')
        pyfits.writeto(fits_file_name, self._intMat, overwrite=True)
        fits_file_name = os.path.join(dove, 'Reconstructor.fits')
        pyfits.writeto(fits_file_name, self._rec, overwrite=True)

    def _saveMask(self, dove, mask):
        """
        args:
            dove = path that indicates where to save the mask file
        """
        fits_file_name = os.path.join(dove, 'Mask.fits')
        pyfits.writeto(fits_file_name, mask.astype(int), overwrite=True)
