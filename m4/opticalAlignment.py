'''
@author: cs
'''

from m4.ground.configuration import Configuration
from m4.utils.zernikeOnM4 import ZernikeOnM4
from m4.opticalCalibration import Calibration
from m4.ground import logger
import numpy as np
import pyfits
import os


class Alignment():
    
    def __init__(self, tt):
        self._tt= tt
        self._cal= Calibration()
        self._zOnM4= ZernikeOnM4()
        self._rec= None
        self._intMat= None
        self._mask= None
        
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "Alignment")
        
    
    def ott_align(self):
        self._intMat, self._rec, self._mask= self._loadAlignmentInfo()
        cmd= self._commandGenerator()
        return cmd
        
    
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
        img= np.ma.masked_array(ima[0], mask= np.invert(ima[1].astype(bool)))
        return img
    
    def _commandGenerator(self):
        self._measureOTTPhaseMap()
        img= self._testAlignment_loadMeasureFromFileFits()
        zernikeVector= self._zernikeCoeff(img)
        cmd= -np.dot(self._rec, zernikeVector)
        return cmd
    
    def _measureOTTPhaseMap(self):
        #salva l'interferogramma
        pass
    
    def _zernikeCoeff(self, img):
        coef, mat= self._zOnM4.zernikeFit(img, np.arange(2,11))
        z= np.array([0,1,2,5,6])
        finalCoef= np.zeros(z.shape[0])
        aa= np.arange(finalCoef.shape[0])
        zipped=zip(aa,z)
        for i,j in zipped:
            finalCoef[i]=coef[j]
        return finalCoef
        
    
    
    
        
        
        
        