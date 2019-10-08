'''
@author: cs
'''

from m4.ground.configuration import Configuration
from m4.utils.zernikeOnM4 import ZernikeOnM4
from m4.utils.opticalCalibration import Opt_Calibration
from m4.ground import logger
import numpy as np
import pyfits
import os


class Opt_Alignment():
    
    def __init__(self, tt):
        self._tt= tt
        self._cal= Opt_Calibration()
        self._zOnM4= ZernikeOnM4()
        self._rec= None
        self._intMat= None
        self._mask= None
        
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "Alignment")
        
    
    def opt_align(self):
        logger.log('Calcolo il comando', 'di allineamento', 'per', self._tt)
        self._intMat, self._rec, self._mask= self._loadAlignmentInfo()
        imgf, imgt= self._measureOTTPhaseMap()
        cmdf= self._commandGenerator(imgf)
        cmdt= self._commandGenerator(imgt)
        self.saveCommand(cmdf, 2)
        return cmdf, cmdt
        
    
    def _loadAlignmentInfo(self):
        self._intMat= self._readInfo('InteractionMatrix.fits')
        self._rec= self._readInfo('Reconstructor.fits')
        self._mask= self._readInfo('Mask.fits')
        return self._intMat, self._rec, self._mask
        
        
    def _readInfo(self, fitsName):
        fold= os.path.join(self._cal._storageFolder(), self._tt)
        file= os.path.join(fold, fitsName)
        hduList= pyfits.open(file)
        info= hduList[0].data
        return info
    
    def _testAlignment_loadMeasureFromFileFits(self):
        fold='/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/Allineamento/20191007_134908/img.fits'
        hduList= pyfits.open(fold)
        ima= hduList[0].data
        imgf= np.ma.masked_array(ima[0], mask= np.invert(ima[1].astype(bool)))
        fold='/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/Allineamento/20191007_135037/img.fits'
        hduList= pyfits.open(fold)
        ima= hduList[0].data
        imgt= np.ma.masked_array(ima[0], mask= np.invert(ima[1].astype(bool)))
        return imgf, imgt
    
    def _commandGenerator(self, img):
        image= np.ma.masked_array(img.data, mask= self._mask)
        zernikeVector= self._zernikeCoeff(image)
        cmd= -np.dot(self._rec, zernikeVector)
        return cmd
    
    def _measureOTTPhaseMap(self):
        #salva l'interferogramma
        imgf, imgt= self._testAlignment_loadMeasureFromFileFits()
        return imgf, imgt
    
    def _zernikeCoeff(self, img):
        coef, mat= self._zOnM4.zernikeFit(img, np.arange(2,11))
        z= np.array([0,1,2,5,6])
        finalCoef= np.zeros(z.shape[0])
        aa= np.arange(finalCoef.shape[0])
        zipped=zip(aa,z)
        for i,j in zipped:
            finalCoef[i]=coef[j]
        return finalCoef
        
    def saveCommand(self, cmd, i):
        '''
        arg:
            cmd= vettore contenente il comando da dare ai gradi di libertà
            i= numero identificativo con cui salvare il nome del comando
        '''
        dove= os.path.join(self._storageFolder(), self._tt)
        if not os.path.exists(dove):
            os.makedirs(dove)
            name= 'Command_%03d.fits' %i
            fitsFileName= os.path.join(dove, name)
            pyfits.writeto(fitsFileName, cmd)
        else:
            name= 'Command_%03d.fits' %i
            fitsFileName= os.path.join(dove, name)
            pyfits.writeto(fitsFileName, cmd)
    
        
        
        
        