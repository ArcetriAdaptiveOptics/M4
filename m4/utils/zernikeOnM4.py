'''
@author: cs
'''

from m4.ground.configuration import Configuration
from m4.ground.zernikeGenerator import ZernikeGenerator
from m4.ground.zernikeMask import CircularMask
import numpy as np


class ZernikeOnM4(object):
    
    def __init__(self):
        self._pupilXYRadius= Configuration.ParabolaPupilXYRadius
        self._zg= ZernikeGenerator(2*self._pupilXYRadius[2])
        
    def getPupilCenterAndRadiusInIFCoords(self):
        return self._pupilXYRadius
        
    def setPupilCenterAndRadiusInIFCoords(self, centerX, centerY, radius):
        self._pupilXYRadius= np.array([centerX, centerY, radius])
        self._zg= ZernikeGenerator(2*radius)
        
    
    def zernikeFit(self, img, zernikeMode):
        '''
        zernikeMode= vector of Zernike modes to remove
        '''
        mat= np.zeros((img.compressed().shape[0], zernikeMode.size)) 
        for i in range(0, zernikeMode.size):
            z=self._zg.getZernike(zernikeMode[i])
            aa= np.ma.masked_array(z, mask= img.mask)
            mat.T[i]= aa.compressed()
            
        self._mat= mat
        inv= np.linalg.pinv(mat)   
        a= np.dot(inv, img.compressed())  
        return a, mat
    
    def zernikeSurface(self, surfaceZernikeCoeffArray, imaMask, mat, index= None):
        zernikeSurfaceMap= None
        
        if index is None:
            zernikeSurfaceMap= np.dot(mat, surfaceZernikeCoeffArray)
        else:
            for i in range(len(index)):
                k= index[i]
                zernikeSurface= np.dot(mat[:,k], surfaceZernikeCoeffArray[k])
                if zernikeSurfaceMap is None:
                    zernikeSurfaceMap= zernikeSurface
                else:
                    zernikeSurfaceMap= zernikeSurfaceMap + zernikeSurface
                
        mask= np.invert(imaMask)
        surf= np.ma.masked_array(np.zeros((2*self._pupilXYRadius[2],2*self._pupilXYRadius[2])), mask=mask)               
        surf[mask]= zernikeSurfaceMap 
        return surf
    
    
    def zernikeToDMCommand(self, surfaceMap, an):
        self._an= an
        '''
        @params
        surfaceMap: 
        an: analyzer delle funzioni d'influenza zonali
        
        @return: 
        Command (numpy.array)
        '''
        cmask = self._an.getMasterMask()
        zernikeMaskOnIF= CircularMask(
            self._an.getIFShape(),
            self._pupilXYRadius[2],
            [self._pupilXYRadius[1], self._pupilXYRadius[0]]).zernikeMask()
        self._an.setDetectorMask(zernikeMaskOnIF | cmask)
        rec=self._an.getReconstructor()
        
        zernikeCmd= np.dot(rec, surfaceMap.compressed())

        return zernikeCmd
    