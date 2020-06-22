'''
@author: cs
'''

import os
import logging
from astropy.io import fits as pyfits
import numpy as np
from m4.ground import tracking_number_folder
from m4.configuration.config import path_name
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from m4.utils.optical_calibration import opt_calibration
from m4.ground import object_from_fits_file_name as obj
from m4.utils.interface_4D import comm4d
from m4.configuration.ott_parameters import OttParameters


class opt_alignment():
    """
    Class for the optical alignment

    HOW TO USE IT::

        from m4.utils.optical_alignement import opt_alignement
        al = opt_alignement()
    """

    def __init__(self, tt):
        """The constructor """
        self._logger = logging.getLogger('OPT_ALIGN:')
        self._tt = tt
        self._cal = opt_calibration()
        self._zOnM4 = ZernikeOnM4()
        self._c4d = comm4d()
        self._rec = None
        self._intMat = None
        self._mask = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save data"""
        return os.path.join(path_name.OPD_DATA_FOLDER,
                            "Alignment")


    def opt_align(self, ott, piston=None):
        """
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
        self._logger.info('Calculation of the alignment command for %s',
                          self._tt)
        self._intMat, self._rec, self._mask = self._loadAlignmentInfo()
        img = self._measureOTTPhaseMap(ott)
        cmd = self._commandGenerator(img)
#         cmdf= self._commandGenerator(imgf)
#         cmdt= self._commandGenerator(imgt)
#         self.saveCommand(cmdf, 1)
        par_command, rm_command = self._reorgCmd(cmd)
        self._saveAllData(par_position, rm_position, par_command, rm_command)
        return par_command, rm_command


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

    def _reorgCmd(self, cmd):
        dofIndex = np.append(OttParameters.PARABOLA_DOF, OttParameters.RM_DOF)
        par_command = np.zeros(6)
        rm_command = np.zeros(6)
        for i in range(cmd.size):
            if i <3:
                par_command[dofIndex[i]] = cmd[i]
            else:
                rm_command[dofIndex[i]] = cmd[i]
        return par_command, rm_command

    def _testAlignment_loadMeasureFromFileFits(self, test):
        """ Test function """
        if test == 0:
            # my pc
            imgf = obj.readImageFromFitsFileName('ProvaCodice/Immagini_prova/Allineamento/20191007_134908/img.fits')
            imgt = obj.readImageFromFitsFileName('ProvaCodice/Immagini_prova/Allineamento/20191007_135037/img.fits')
            # m4 pc
#            imgf = obj.readImageFromFitsFileName('TestData/Allineamento/20191007_134908/img.fits')
#            imgt = obj.readImageFromFitsFileName('TestData/Allineamento/20191007_135037/img.fits')
            return imgf, imgt
        elif test == 1:
            # my pc
            img = obj.readImageFromFitsFileName('ProvaCodice/Immagini_prova/Allineamento/20191001_081344/img.fits')
            # m4 pc
#            img = obj.readImageFromFitsFileName('TestData/Allineamento/20191001_081344/img.fits')
            return img

    def _testAlignment_loadInfoPositionFromFileFits(self, test):
        """ Test function """
        if test == 0:
            # my pc
#             parf = obj.readDataFromFileFits('Immagini_prova/Allineamento/20191007_134908/parpos.fits')
#             rmf = obj.readDataFromFileFits('Immagini_prova/Allineamento/20191007_134908/rflatpos.fits')
#             part = obj.readDataFromFileFits('Immagini_prova/Allineamento/20191007_135037/parpos.fits')
#             rmt = obj.readDataFromFileFits('Immagini_prova/Allineamento/20191007_135037/rflatpos.fits')
            # m4 pc
            parf = obj.readDataFromFileFits('TestData/Allineamento/20191007_134908/parpos.fits')
            rmf = obj.readDataFromFileFits('TestData/Allineamento/20191007_134908/rflatpos.fits')
            part = obj.readDataFromFileFits('TestData/Allineamento/20191007_135037/parpos.fits')
            rmt = obj.readDataFromFileFits('TestData/Allineamento/20191007_135037/rflatpos.fits')
            return parf, rmf, part, rmt
        elif test == 1:
            # my pc
            par = obj.readDataFromFileFits('Immagini_prova/Allineamento/20191001_081344/parpos.fits')
            rm = obj.readDataFromFileFits('Immagini_prova/Allineamento/20191001_081344/rflatpos.fits')
            # m4 pc
#            par = obj.readDataFromFileFits('TestData/Allineamento/20191001_081344/parpos.fits')
#            rm = obj.readDataFromFileFits('TestData/Allineamento/20191001_081344/rflatpos.fits')
            return par, rm

    def _commandGenerator(self, img):
        """
        args:
            img = image

        returns:
                cmd = command for the dof
        """
        image = np.ma.masked_array(img.data, mask=self._mask)
        zernike_vector = self._zernikeCoeff(image)
        cmd = - np.dot(self._rec, zernike_vector)
        return cmd

    def _measureOTTPhaseMap(self, ott):
        #acquisirà e salverà l'interferogramma
        self._logger.debug('Measure of phase map')
#         imgf, imgt = self._testAlignment_loadMeasureFromFileFits(0)
#         img = self._testAlignment_loadMeasureFromFileFits(1)
        p, m = self._c4d.acq4d(ott, 1, show=1)
        img = np.ma.masked_array(p, mask=np.invert(m.astype(bool)))
        return img

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

    def _saveAllData(self, par_position, rm_position, par_command, rm_command):
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


    def saveCommand(self, cmd, i):
        '''
        Parameters
        ----------
            cmd : numpy array
                vector containing the command to be given to degrees of freedom
            i: int
                id number to save the command name
        '''
        dove = os.path.join(self._storageFolder(), self._tt)
        if not os.path.exists(dove):
            os.makedirs(dove)
            name = 'Command_%03d.fits' %i
            fits_file_name = os.path.join(dove, name)
            pyfits.writeto(fits_file_name, cmd)
        else:
            name = 'Command_%03d.fits' %i
            fits_file_name = os.path.join(dove, name)
            pyfits.writeto(fits_file_name, cmd)
