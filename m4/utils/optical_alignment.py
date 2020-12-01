'''
Autors
  - C. Selmi:  written in 2019
'''

import os
import logging
from astropy.io import fits as pyfits
import numpy as np
from m4.ground import tracking_number_folder
from m4.configuration.config import fold_name
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from m4.utils.optical_calibration import opt_calibration
from m4.utils.interface_4D import comm4d
from m4.configuration.ott_parameters import OttParameters
from m4.configuration import config as conf
from m4.utils import image_extender as ie


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
        self._zOnM4 = ZernikeOnM4()
        self._c4d = comm4d()
        self._rec = None
        self._intMat = None
        self._mask = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save data"""
        return fold_name.ALIGNMENT_ROOT_FOLDER


    def opt_align(self, ott, n_images, piston=None):
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
        self._intMat, self._rec, self._mask = self._loadAlignmentInfo()
        img = self._measureOTTPhaseMap(ott, n_images)

        if self._cal._who=='PAR + RM':
            cmd = self._commandGenerator(img)
            par_command, rm_command = self._reorgCmdMix(cmd)
            self._saveAllDataMix(par_position, rm_position, par_command, rm_command)
            return par_command, rm_command
        elif self._cal._who=='M4':
            cmd = self._commandGenerator(img, piston)
            m4_command = self._reorgCmdM4(cmd)
            self._saveAllDataM4(m4_position, m4_command)
            return m4_command


    def _loadAlignmentInfo(self):
        """ Returns interaction matrix, reconstructor and mask """
        self._intMat = self._readInfo('InteractionMatrix.fits')
        self._rec = self._readInfo('Reconstructor.fits')
        self._mask = self._readInfo('Mask.fits')
        return self._intMat, self._rec, self._mask


    def _readInfo(self, fits_name):
        """ Function for reading fits file"""
        fold = os.path.join(self._cal._storageFolder(), self._tt)
        file = os.path.join(fold, fits_name)
        hduList = pyfits.open(file)
        info = hduList[0].data
        return info

    def _reorgCmdMix(self, cmd):
        dofIndex = np.append(OttParameters.PARABOLA_DOF, OttParameters.RM_DOF)
        par_command = np.zeros(6)
        rm_command = np.zeros(6)
        for i in range(cmd.size):
            if i < OttParameters.PARABOLA_DOF.size:
                par_command[dofIndex[i]] = cmd[i]
            else:
                rm_command[dofIndex[i]] = cmd[i]
        return par_command, rm_command

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
        #image = np.ma.masked_array(img.data, mask=self._mask)
        new_image = ie.imageExtender(img)
        zernike_vector = self._zernikeCoeff(new_image)
        print('zernike:')
        print(zernike_vector)
        if piston is None:
            zernike_vector = zernike_vector
        else:
            cc = zernike_vector[3]
            zernike_vector[3] = cc + piston
        #sommare il coma a questo zernike vector
        cmd = - np.dot(self._rec, zernike_vector)
        print('mix command:')
        print(cmd)
        return cmd

    def _measureOTTPhaseMap(self, ott, n_images):
        #acquisir e salver l'interferogramma
        self._logger.debug('Measure of phase map')
#         imgf, imgt = self._testAlignment_loadMeasureFromFileFits(0)
#         img = self._testAlignment_loadMeasureFromFileFits(1)
        cube_images = None
        for i in range(n_images):
            masked_ima = self._c4d.acq4d(ott, 1, show=0)
            if cube_images is None:
                cube_images = masked_ima
            else:
                cube_images = np.ma.dstack((cube_images, masked_ima))
        final_ima = np.ma.mean(cube_images, axis=2)
        return final_ima

    def _zernikeCoeff(self, img):
        """
        Returns:
                final_coef = zernike coeff on the image
                            (zernike modes 2,3,4,7,8)
        """
        coef, mat = self._zOnM4.zernikeFit(img, np.arange(2, 11))
        z = np.array([0, 1, 2, 5, 6])
        final_coef = np.zeros(z.shape[0])
        aa = np.arange(final_coef.shape[0])
        zipped = zip(aa, z)
        for i, j in zipped:
            final_coef[i] = coef[j]
        return final_coef

    def _saveAllDataMix(self, par_position, rm_position, par_command, rm_command):
        save = tracking_number_folder.TtFolder(self._storageFolder())
        dove, self._align_tt = save._createFolderToStoreMeasurements()
        name = 'par_position.fits'
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, par_position)
        name = 'rm_position.fits'
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, rm_position)
        name = 'par_command.fits'
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, par_command)
        name = 'rm_command.fits'
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, rm_command)

    def _saveAllDataM4(self, m4_position, m4_command):
        save = tracking_number_folder.TtFolder(self._storageFolder())
        dove, self._align_tt = save._createFolderToStoreMeasurements()
        name = 'm4_position.fits'
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, m4_position)
        name = 'm4_command.fits'
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, m4_command)
