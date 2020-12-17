'''
Authors
  - C. Selmi:  written in 2019
'''

import os
import logging
from astropy.io import fits as pyfits
import numpy as np
from m4.configuration.config import fold_name
from m4.utils.optical_calibration import opt_calibration
from m4.utils.interface_4D import comm4d
from m4.ground import zernike
from m4.configuration.ott_parameters import OttParameters
from m4.ground.timestamp import Timestamp
from matplotlib import pyplot as plt


class opt_alignment():
    """
    Class for the optical alignment

    HOW TO USE IT::

        from m4.utils.optical_alignment import opt_alignment
        al = opt_alignment(tt)
        command = al.opt_align(ott, piston=None)
    """

    def __init__(self, tt):
        """The constructor """
        self._logger = logging.getLogger('OPT_ALIGN:')
        self._tt = tt
        self._cal = opt_calibration.loadCommandMatrixFromFits(tt)
        self._c4d = comm4d()
        self._rec = None
        self._intMat = None
        self._mask = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save data"""
        return fold_name.ALIGNMENT_ROOT_FOLDER


    def opt_align(self, ott, n_images, intMatModesVector=None, commandId=None, piston=None):
        """
        Parameters
        ----------
        ott: object
            test tower

        Other Parameters
        ----------
            piston: int, optional

        Returns
        -------
                cmd: numpy array
                     final command for the optical alignment
        """
        par_position = ott.parab()
        rm_position = ott.refflat()
        #m4_position = ott.m4()
        self._logger.info('Calculation of the alignment command for %s',
                          self._tt)
        self._intMat, self._rec, self._mask = self._selectModesInIntMatAndRecConstruction(intMatModesVector, commandId)
        self._intMatModesVector = intMatModesVector

        img = self._c4d.acq4d(ott, n_images, 0)
        name = 'StartImage.fits'
        tt = Timestamp.now()
        dove = os.path.join(self._storageFolder(), self._tt + '--' + tt)
        os.makedirs(dove)
        self._c4d.save_phasemap(dove, name, img)

        if self._cal._who=='PAR + RM':
            cmd, zernike_vector = self._commandGenerator(img)
            par_command, rm_command = self._reorgCmdMix(cmd, commandId)
            self._saveAllDataMix(dove, par_position, rm_position, par_command, rm_command)
            self._saveZernikeVector(dove, zernike_vector)
            return par_command, rm_command, dove
        elif self._cal._who=='M4':
            cmd, zernike_vector = self._commandGenerator(img, piston)
            m4_command = self._reorgCmdM4(cmd)
            self._saveAllDataM4(dove, m4_position, m4_command)
            self._saveZernikeVector(dove, zernike_vector)
            return m4_command, dove

    def _selectModesInIntMatAndRecConstruction(self, intMatModesVector=None,
                                               commandId=None):
        intMat, rec, mask = self._loadAlignmentInfo()
        if intMatModesVector is None:
            new_intMat = intMat
        else:
            new_intMat = intMat[intMatModesVector, :]

        if commandId is not None:
            new_intMat = new_intMat[:,commandId]

        new_rec = np.linalg.pinv(new_intMat)

        return new_intMat, new_rec, mask

    def _loadAlignmentInfo(self):
        """ Returns interaction matrix, reconstructor and mask """
        self._intMat = self._readInfo('InteractionMatrix.fits')
        self._rec = self._readInfo('Reconstructor.fits')
        self._mask = self._readInfo('Mask.fits')
        #y = ['PAR_PIST', 'PAR_TIP', 'PAR_TILT', 'RM_TIP', 'RM_TILT']
        plt.clf()
        plt.imshow(self._intMat, origin='lower')
        plt.colorbar()
        plt.xlabel('Commands')
        plt.ylabel('Zernike Modes')
        return self._intMat, self._rec, self._mask


    def _readInfo(self, fits_name):
        """ Function for reading fits file"""
        fold = os.path.join(self._cal._storageFolder(), self._tt)
        file = os.path.join(fold, fits_name)
        hduList = pyfits.open(file)
        info = hduList[0].data
        return info

    def _reorgCmdMix(self, cmd, commandId=None):
        dofIndex = np.append(OttParameters.PARABOLA_DOF, OttParameters.RM_DOF)
        par_command = np.zeros(6)
        rm_command = np.zeros(6)

        if commandId is not None:
            mycomm = np.zeros(5)
            mycomm[commandId]=cmd
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

    def _commandGenerator(self, img, piston=None):
        """
        args:
            img = image

        returns:
                cmd = command for the dof
        """
        zernike_vector = self._zernikeCoeff(img)
        print('zernike:')
        print(zernike_vector)
        if piston is None:
            zernike_vector = zernike_vector
        else:
            cc = zernike_vector[3]
            zernike_vector[3] = cc + piston
        #sommare il coma a questo zernike vector
        cmd = - np.dot(self._rec, zernike_vector) #giusto
        print('mix command:')
        print(cmd)
        return cmd, zernike_vector


    def _zernikeCoeff(self, img):
        """
        Returns:
                final_coef = zernike coeff on the image
                            (zernike modes 2,3,4,7,8)
        """
        coef, mat = self._zOnM4.zernikeFit(img, np.arange(10)+1)
        z = np.array([1, 2, 3, 6, 7])
        final_coef = coef[z]

        if self._intMatModesVector is None:
            final_coef_selected = final_coef
        else:
            final_coef_selected = np.zeros(self._intMatModesVector.size)
            for i in range(self._intMatModesVector.size):
                final_coef_selected[i] = final_coef[self._intMatModesVector[i]]
        return final_coef_selected

    def _saveAllDataMix(self, dove, par_position, rm_position, par_command, rm_command):
        name = 'PositionAndDeltaCommand.fits'
        vector = np.array([par_position, rm_position, par_command, rm_command])
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, vector)

    def _saveAllDataM4(self, dove, m4_position, m4_command):
        name = 'PositionAndDeltaCommand.fits'
        vector = np.array([m4_position, m4_command])
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, vector)

    def _saveZernikeVector(self, dove, zernike_vector):
        name = 'Zernike.fits'
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, zernike_vector)
