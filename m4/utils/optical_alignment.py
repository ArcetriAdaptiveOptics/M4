'''
Authors
  - C. Selmi:  written in 2019
'''

import os
import logging
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.configuration.config import fold_name
from m4.configuration import config as conf
from m4.utils.optical_calibration import OpticalCalibration
from m4.ground import zernike
from m4.configuration.ott_parameters import OttParameters, OtherParameters
from m4.ground.timestamp import Timestamp
from m4.utils.roi import ROI


class OpticalAlignment():
    """
    Class for the optical alignment

    HOW TO USE IT::

        from m4.utils.optical_alignment import opt_alignment
        al = OpticalAlignment(tt)
        command = al.opt_align(ott, piston=None)
    """

    def __init__(self, tt, ott, interf):
        """The constructor """
        self._logger = logging.getLogger('OPT_ALIGN:')
        self._tt = tt
        self._cal = OpticalCalibration.loadCommandMatrixFromFits(tt)
        self._interf = interf
        self._ott = ott
        self._rec = None
        self._intMat = None
        self._mask = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save data"""
        return fold_name.ALIGNMENT_ROOT_FOLDER


    def opt_align(self, n_images, intMatModesVector=None, commandId=None):
        """
        Parameters
        ----------
        ott: object
            test tower
        n_images: int
            number of interferometers frames

        Other Parameters
        ----------
            intMatModesVecor: numpy array
                        None is equal to np.array([0,1,2,3,4,5])
                        for tip, tilt, fuoco, coma, coma
            commandId: numpy array
                    array containing the number of degrees of freedom to be commanded
            piston: int, optional

        Returns
        -------
                cmd: numpy array
                    final delta command for the optical alignment
        """
        par_position = self._ott.parabola.getPosition()
        rm_position = self._ott.referenceMirror.getPosition()
        m4_position = self._ott.m4.getPosition()
        self._logger.info('Calculation of the alignment command for %s',
                          self._tt)
        self._intMat, self._rec, self._cmat, self._mask = self._selectModesInIntMatAndRecConstruction(intMatModesVector, commandId)
        self._intMatModesVector = intMatModesVector

        img = self._interf.acquire_phasemap(n_images)
        name = 'StartImage.fits'
        tt = Timestamp.now()
        dove = os.path.join(self._storageFolder(), self._tt + '--' + tt)
        os.makedirs(dove)
        self._interf.save_phasemap(dove, name, img)

        if self._cal._who=='PAR + RM':
            cmd, zernike_vector, total_coef = self._commandGenerator(img)
            par_command, rm_command = self._reorgCmdMix(cmd, commandId)
            self._saveAllDataMix(dove, par_position, rm_position, par_command, rm_command,
                                 intMatModesVector, commandId)
            self._saveZernikeVector(dove, zernike_vector)
            self._alignmentLog(total_coef, tt, commandId)
            return par_command, rm_command, dove
        elif self._cal._who=='M4':
            cmd, zernike_vector = self._commandGenerator(img)
            m4_command = self._reorgCmdM4(cmd)
            self._saveAllDataM4(dove, m4_position, m4_command)
            self._saveZernikeVector(dove, zernike_vector)
            return m4_command, dove

    def _alignmentLog(self, start_total_coef, tt, commandId):
        fits_file_name = os.path.join(self._storageFolder(), 'AlignmentLog.txt')
        file = open(fits_file_name, 'a+')
        file.write('%s ' %self._tt)
        for i in range(start_total_coef.size):
            file.write('%9.3e ' %start_total_coef[i])
        file.write('\n')
        file.write('%s ' %tt)
#         for i in range(total_coef.size):
#             file.write('%9.3e ' %total_coef[i])
#         file.write('\n')
#         file.write('%s \n ************\n' %commandId)
        file.close()

    def _selectModesInIntMatAndRecConstruction(self, zernike2control=None,
                                               commandId=None):
        '''
        Other Parameters
        ----------
            intMatModesVecor: numpy array
                        None is equal to np.array([0,1,2,3,4,5])
                        for tip, tilt, fuoco, coma, coma
            commandId: numpy array
                    array containing the number of degrees of freedom to be commanded
        '''
        intMat, mask = self._loadAlignmentInfo()
        cmat = self._cal._commandMatrix
        if zernike2control is None:
            new_intMat = intMat
            new_cmat = cmat
        else:
            new_intMat = intMat[zernike2control, :]

        if commandId is not None:
            new_intMat = new_intMat[:, commandId]
            new_cmat = cmat[commandId, :]
            new_cmat = new_cmat[:, commandId]

        new_rec = np.linalg.pinv(new_intMat)
        return new_intMat, new_rec, new_cmat, mask

    def _loadAlignmentInfo(self):
        """ Returns interaction matrix, reconstructor and mask """
        self._intMat = self._readInfo('InteractionMatrix.fits')
        self._mask = self._readInfo('Mask.fits')
        #y = ['PAR_PIST', 'PAR_TIP', 'PAR_TILT', 'RM_TIP', 'RM_TILT']
        plt.clf()
        plt.imshow(self._intMat, origin='lower')
        plt.colorbar()
        plt.xlabel('Commands')
        plt.ylabel('Zernike Modes')
        return self._intMat, self._mask


    def _readInfo(self, fits_name):
        """ Function for reading fits file"""
        fold = os.path.join(self._cal._storageFolder(), self._tt)
        file = os.path.join(fold, fits_name)
        hduList = pyfits.open(file)
        info = hduList[0].data
        return info

    def _reorgCmdMix(self, cmd, commandId=None):
        '''reorganizes the delta command in the
        right positions for par and rm '''
        dofIndex = np.append(OttParameters.PARABOLA_DOF, OttParameters.RM_DOF)
        par_command = np.zeros(6)
        rm_command = np.zeros(6)

        if commandId is not None:
            mycomm = np.zeros(5)
            mycomm[commandId] = cmd
            cmd = mycomm

        for i in range(cmd.size):
            if i < OttParameters.PARABOLA_DOF.size:
                par_command[dofIndex[i]] = cmd[i]
            else:
                rm_command[dofIndex[i]] = cmd[i]

        return par_command, rm_command

#va riscritta
    def _reorgCmdM4(self, cmd):
        dofIndex = OttParameters.M4_DOF
        m4_command = np.zeros(6)
        for i in range(cmd.size):
            m4_command[dofIndex[i]] = cmd[i]
        return m4_command

    def _commandGenerator(self, img):
        """
        args:
            img = image

        returns:
                cmd = command for the dof
        """
        total_coef, zernike_vector = self._zernikeCoeff(img)
        print('zernike:')
        print(zernike_vector)
#         if old_or_new==1:
#             cmd = - np.dot(self._rec, zernike_vector)
#         else:
        M = np.dot(self._cmat, self._rec) #non serve la trasposta
        cmd = - np.dot(M, zernike_vector) #serve il meno
        print('mix command:')
        print(cmd)
        return cmd, zernike_vector, total_coef


    def _zernikeCoeff(self, img):
        """
        Returns:
                final_coef = zernike coeff on the image
                            (zernike modes 2,3,4,7,8)
        """
        if  conf.simulated ==1:
            mask_index = OtherParameters.MASK_INDEX_SIMULATORE
        else:
            mask_index = OtherParameters.MASK_INDEX_TOWER
        r = ROI()
        roi = r.roiGenerator(img)
        mask = roi[mask_index]
        mm = np.ma.mask_or(img.mask, mask)

        new_image = np.ma.masked_array(img, mask=mm)
        coef, mat = zernike.zernikeFit(new_image, np.arange(10)+1)
        z = np.array([1, 2, 3, 6, 7])
        final_coef = coef[z]

        if self._intMatModesVector is None:
            final_coef_selected = final_coef
        else:
            final_coef_selected = np.zeros(self._intMatModesVector.size)
            for i in range(self._intMatModesVector.size):
                final_coef_selected[i] = final_coef[self._intMatModesVector[i]]
        return final_coef, final_coef_selected

    def _saveAllDataMix(self, dove, par_position, rm_position,
                        par_command, rm_command, intMatModesVector,
                        commandId):
        name = 'PositionAndDeltaCommand.fits'
        vector = np.array([par_position, rm_position, par_command, rm_command])
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, vector)
        if intMatModesVector is not None:
            fits_file_name = os.path.join(dove, 'intMatModesVector')
            pyfits.writeto(fits_file_name, intMatModesVector)
        if commandId is not None:
            fits_file_name = os.path.join(dove, 'commandId')
            pyfits.writeto(fits_file_name, commandId)

    def _saveAllDataM4(self, dove, m4_position, m4_command):
        name = 'PositionAndDeltaCommand.fits'
        vector = np.array([m4_position, m4_command])
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, vector)

    def _saveZernikeVector(self, dove, zernike_vector):
        name = 'Zernike.fits'
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, zernike_vector)
