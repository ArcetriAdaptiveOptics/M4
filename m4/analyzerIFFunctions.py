'''
@author: cs
'''

from m4.ground import logger
import numpy as np
import pyfits
import os
from m4.ground.interferometer_converter import InterferometerConverter
from m4.influenceFunctionsMaker import IFFunctionsMaker

class AnalyzerIFF(object):
    
    def __init__(self):
        self._indexingList= None
        self._cube= None
        self._rec= None
        self._intMat= None
        self._analysisMask= None

    def _ttData(self):
        split= os.path.split(self._h5Folder)
        self._tt= split[1]
        return self._tt
    
    @staticmethod
    def loadInfoFromh5Folder(h5Folder):
        theObject= AnalyzerIFF()
        theObject._h5Folder= h5Folder
        a= IFFunctionsMaker.loadInfo(h5Folder)
        theObject._who= a[0]
        theObject._actsVector= a[1]
        theObject._cmdMatrix= a[2]
        theObject._cmdAmplitude= a[3]
        theObject._nPushPull= a[4]
        theObject._indexingList= a[5]
        return theObject
    
    
    
    def getCube(self):
        return self._cube
    
    def getIFShape(self):
        return self.getCube()[:,:,0].shape
    
    def getMasterMask(self):
        aa=np.sum(self._cube.mask.astype(int), axis=2)
        masterMask= np.zeros(aa.shape, dtype=np.bool)
        masterMask[np.where(aa>0)]= True
        return masterMask    
    
    
    def _indexReorganization(self):
        indv= np.array(self._indexingList)
        where=[]
        for ind in self._actsVector:
            for j in range(self._nPushPull):
                a=np.where(indv[j]==ind) 
                where.append(a)      
        return where
    
    
    def createCube(self):
        cubeAllAct= None
        ic= InterferometerConverter()
        self._ttData()
        logger.log('Creazione del cubo delle IFF per', self._who, self._tt)
        

        self._cube= cubeAllAct                
        return 
    
    
    def saveCubeAsFits(self, fitsFileName):
        header= pyfits.Header()
        header['NPUSHPUL']= self._nPushPull
        header['WHO']= self._who
        pyfits.writeto(fitsFileName, self._cube.data, header)
        pyfits.append(fitsFileName, self._cube.mask.astype(int))
        pyfits.append(fitsFileName, self._cmdAmplitude)
        pyfits.append(fitsFileName, self._actsVector)
        
    @staticmethod
    def loadCubeFromFits(fitsFileName):
        header= pyfits.getheader(fitsFileName)
        hduList= pyfits.open(fitsFileName)
        cube= np.ma.masked_array(hduList[0].data, hduList[1].data.astype(bool))
        actsVector= hduList[3].data
        cmdAmplitude= hduList[2].data
        who= header['WHO']
        try:
            nPushPull= header['NPUSHPUL']
        except KeyError:
            nPushPull= 1
        
        theObject= AnalyzerIFF()
        theObject._actsVector= actsVector
        theObject._cmdAmplitude= cmdAmplitude
        theObject._nPushPull= nPushPull
        theObject._who= who
        theObject._cube= cube
        return theObject
    
    
    def setAnalysisMask(self, analysisMask):
        self._analysisMask= analysisMask
    
    def setAnalysisMaskFromMasterMask(self):
        self._analysisMask= self.getMasterMask()
        
    def setDetectorMask(self, maskFromIma):
        self._analysisMask= maskFromIma
        self._rec=None
        
    def _getMaskedInfluenceFunction(self, idxInfluenceFunction):
        return np.ma.array(self.getCube()[:,:,idxInfluenceFunction], 
                            mask=self.getAnalysisMask())
    
    def _createInteractionMatrix(self):
        if self._analysisMask is None:
            self.setAnalysisMaskFromMasterMask()
        nActsInCube= self.getCube().shape[2] 
        nInterferometerPixelsInMask= self._getMaskedInfluenceFunction(0).compressed().shape[0]
        self._intMat=np.zeros((nInterferometerPixelsInMask, nActsInCube))
        for i in range(nActsInCube):                                          
            self._intMat[:, i]= self._getMaskedInfluenceFunction(i).compressed()
    
    def _createSurfaceReconstructor(self, rCond= 1e-15):
        self._rec= self._createRecWithPseudoInverse(rCond)
    
    def _createRecWithPseudoInverse(self, rCond):
        return np.linalg.pinv(self.getInteractionMatrix(), rcond=rCond)
    
    def getInteractionMatrix(self):
        if self._intMat is None:
            self._createInteractionMatrix()
        return self._intMat
    
    def getReconstructor(self):
        if self._rec is None:
            self._createSurfaceReconstructor()
        return self._rec
    
    
    
    