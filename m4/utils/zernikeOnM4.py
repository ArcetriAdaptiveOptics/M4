'''
@author: cs
'''

from m4.ground.zernikeGenerator import ZernikeGenerator
from m4.ground.zernikeMask import CircularMask
import numpy as np


class ZernikeOnM4(object):
    
    def __init__(self, analyzerIFFunctions):
        self._an = analyzerIFFunctions
        self._shapeIFs= self._an.getCube()[:,:,0].shape
        self._nPixelOnDiameter= None
        self._nActs= self._an.getCube()[:,0,0].shape
        self._pupilXYRadius= None
        self._zg= None
        
        
    def setPupilCenterAndRadiusInIFCoords(self, centerX, centerY, radius):
        self._pupilXYRadius= np.array([centerX, centerY, radius])
        self._zg= ZernikeGenerator(2*radius)
        
    
    def zernikeFit(self, img, zernikeMode):
        '''
        zernikeMode= vector of Zernike modes to remove
        '''
        z= self._zg.getZernike(2)
        nPointsInMask= len(z.compressed())
        cx= self._pupilXYRadius[0]
        cy= self._pupilXYRadius[1]
        r= self._pupilXYRadius[2]
        imaCut=img[cy-r:cy+r, cx-r:cx+r]
        imaCutM= np.ma.masked_array(imaCut.data, mask=z.mask)
        
        a=np.zeros(zernikeMode.size)
        
        for i in range(0, zernikeMode.size):
            z=self._zg.getZernike(zernikeMode[i])
            a[i]= np.dot(z.compressed(),imaCutM.compressed()) / nPointsInMask
            
        return a
    
    
    def zernikeToDMCommand(self, surfaceZernikeCoeffArrayInNanometer):
        '''
        @params
        surfaceZernikeCoeffArrayInNanometer: numpy.array() coefficients in meter
          Starting from z2
        
        @return: 
        Command (numpy.array)
        '''
        surfaceMap=0.0;
        firstZernModeIndex= 2
        lastZernModeIndex= 2+len(surfaceZernikeCoeffArrayInNanometer)
        indexZernModes= np.arange(firstZernModeIndex, lastZernModeIndex)
        zd= self._zg.getZernikeDict(indexZernModes)
            
        for i in indexZernModes:
            surfaceMap= surfaceMap+ surfaceZernikeCoeffArrayInNanometer[i-2]*zd[i]
         
         
        
        cmask = self._an.getMasterMask()
        zernikeMaskOnIF= CircularMask(
            self._an.getIFShape(),
            self._pupilXYRadius[2],
            [self._pupilXYRadius[1], self._pupilXYRadius[0]]).zernikeMask()
        self._an.setDetectorMask(zernikeMaskOnIF | cmask)
        rec=self._an.getReconstructor()
        
        zernikeCmd= np.dot(rec, surfaceMap.compressed())

        return zernikeCmd, surfaceMap
    