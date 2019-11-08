'''
@author: cs
'''

import os
import logging
import pyfits
import numpy as np
from m4.ground.configuration import Configuration
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from m4.utils.optical_calibration import opt_calibration
from m4.ground import object_from_fits_file_name as obj


class opt_alignment():

    def __init__(self, tt):
        self._logger = logging.getLogger('OPT_ALIGN:')
        self._tt = tt
        self._cal = opt_calibration()
        self._zOnM4 = ZernikeOnM4()
        self._rec = None
        self._intMat = None
        self._mask = None

    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "Alignment")


    def opt_align(self, piston=None):
        self._logger.info('Calculation of the alignment command for %s', self._tt)
        self._intMat, self._rec, self._mask = self._loadAlignmentInfo()
        img = self._measureOTTPhaseMap()
        cmd = self._commandGenerator(img)
#         cmdf= self._commandGenerator(imgf)
#         cmdt= self._commandGenerator(imgt)
#         self.saveCommand(cmdf, 1)
        return cmd


    def _loadAlignmentInfo(self):
        self._intMat = self._readInfo('InteractionMatrix.fits')
        self._rec = self._readInfo('Reconstructor.fits')
        self._mask = self._readInfo('Mask.fits')
        return self._intMat, self._rec, self._mask


    def _readInfo(self, fits_name):
        fold = os.path.join(self._cal._storageFolder(), self._tt)
        file = os.path.join(fold, fits_name)
        hduList = pyfits.open(file)
        info = hduList[0].data
        return info

    def _testAlignment_loadMeasureFromFileFits(self, test):
        if test == 0:
            imgf = obj.readImageFromFitsFileName('Allineamento/20191007_134908/img.fits')
            imgt = obj.readImageFromFitsFileName('Allineamento/20191007_135037/img.fits')
            return imgf, imgt
        elif test == 1:
            img = obj.readImageFromFitsFileName('Allineamento/20191001_081344/img.fits')
            return img

    def _testAlignment_loadInfoPositionFromFileFits(self, test):
        if test == 0:
            parf = obj.readDataFromFileFits('Allineamento/20191007_134908/parpos.fits')
            rmf = obj.readDataFromFileFits('Allineamento/20191007_134908/rflatpos.fits')
            part = obj.readDataFromFileFits('Allineamento/20191007_135037/parpos.fits')
            rmt = obj.readDataFromFileFits('Allineamento/20191007_135037/rflatpos.fits')
            return parf, rmf, part, rmt
        elif test == 1:
            par = obj.readDataFromFileFits('Allineamento/20191001_081344/parpos.fits')
            rm = obj.readDataFromFileFits('Allineamento/20191001_081344/rflatpos.fits')
            return par, rm

    def _commandGenerator(self, img):
        image = np.ma.masked_array(img.data, mask=self._mask)
        zernike_vector = self._zernikeCoeff(image)
        cmd = - np.dot(self._rec, zernike_vector)
        return cmd

    def _measureOTTPhaseMap(self):
        #acquisirà e salverà l'interferogramma
        self._logger.debug('Measure of phase map')
        imgf, imgt = self._testAlignment_loadMeasureFromFileFits(0)
        img = self._testAlignment_loadMeasureFromFileFits(1)
        return img

    def _zernikeCoeff(self, img):
        coef, mat = self._zOnM4.zernikeFit(img, np.arange(2, 11))
        z = np.array([0, 1, 2, 5, 6])
        final_coef = np.zeros(z.shape[0])
        aa = np.arange(final_coef.shape[0])
        zipped = zip(aa, z)
        for i, j in zipped:
            final_coef[i] = coef[j]
        return final_coef

    def saveCommand(self, cmd, i):
        '''
        arg:
            cmd= vettore contenente il comando da dare ai gradi di libertà
            i= numero identificativo con cui salvare il nome del comando
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
