'''
@author: cs
'''

from m4.ground import logger
import numpy as np
import pyfits
import os
from m4.ground.interferometer_converter import InterferometerConverter
from m4.influenceFunctionsMaker import IFFunctionsMaker
from pylint.test.functional import bad_staticmethod_argument
from _pylief import NONE

class AnalyzerIFF(object):
    
    def __init__(self):
        self._indexingList= None
        self._cube= None
        self._rec= None
        self._intMat= None
        self._analysisMask= None
        self._cubeMeasure= None

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
    
    
    def createCube(self, tt):
        cubeAllAct= None
        #self._ttData()
        #logger.log('Creazione del cubo delle IFF per', self._who, self._tt)
        self._cubeMeasure, self._actsVector, cmdMatrix, cmdAmplitude, self._indexingList, who, tt_cmdH, self._nPushPull= self.loadMeasureFromFits(tt)
        where= self._indexReorganization()
        misurePos, misureNeg= self._splitMeasureFromFits(self._cubeMeasure)
        
        for i in range(self._actsVector.shape[0]):
            for k in range(self._nPushPull):
                p= self._nPushPull * i + k
                where= self._indexReorganization()
                n=where[p][0][0]
                mis= k* self._indexingList.shape[1] + n
                
                imgPos= misurePos[:,:,mis]
                imgNeg= misureNeg[:,:,mis]
                imgIF= (imgPos - imgNeg) / 2 #non divido per l'ampiezza perchè nelle misure non c'è
                ifPushPullKth= imgIF-np.ma.median(imgIF)
                
                if k==0:
                    allPushPullActJth= ifPushPullKth
                else:
                    allPushPullActJth= np.ma.dstack((allPushPullActJth, ifPushPullKth))
                    
            if self._nPushPull==1:
                ifActJth= allPushPullActJth
            else:
                ifActJth= np.ma.mean(allPushPullActJth, axis=2)   
                
            if cubeAllAct is None:
                cubeAllAct= ifActJth
            else:
                cubeAllAct= np.ma.dstack((cubeAllAct, ifActJth))
        self._cube= cubeAllAct
                               
        return self._cube
    
        
    def _splitMeasureFromFits(self, misure):
        misurePos= None 
        misureNeg= None 
        for j in range(misure.shape[2]):
            if j%2 == 0:
                if misurePos is None:
                    misurePos= misure[:,:,j]
                else:
                    misurePos= np.ma.dstack((misurePos, misure[:,:,j]))
            else:
                if misureNeg is None:
                    misureNeg= misure[:,:,j]
                else:
                    misureNeg= np.ma.dstack((misureNeg, misure[:,:,j]))
        
        return misurePos, misureNeg
    
    
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
    
    @staticmethod 
    def loadMeasureFromFits(tt):
        dove= os.path.join("/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions", tt)
        fitsFileName= os.path.join(dove, 'misure.fits')
        header= pyfits.getheader(fitsFileName)
        hduList= pyfits.open(fitsFileName)
        misure= np.ma.masked_array(hduList[4].data, hduList[5].data.astype(bool))
        actsVector= hduList[0].data
        cmdMatrix= hduList[1].data
        cmdAmplitude= hduList[2].data
        indexingList= hduList[3].data
        who= header['WHO']
        tt_cmdH= header['TT_CMDH']
        try:
            nPushPull= header['NPUSHPUL']
        except KeyError:
            nPushPull= 1
        
        return misure, actsVector, cmdMatrix, cmdAmplitude, indexingList, who, tt_cmdH, nPushPull

        
    
    def setAnalysisMask(self, analysisMask):
        self._analysisMask= analysisMask
    
    def setAnalysisMaskFromMasterMask(self):
        self._analysisMask= self.getMasterMask()
        
    def setDetectorMask(self, maskFromIma):
        self._analysisMask= maskFromIma
        self._rec=None
        
    def getAnalysisMask(self):
        return self._analysisMask
        
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
    
    
    
    