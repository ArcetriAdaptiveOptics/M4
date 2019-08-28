'''
@author: cs
'''
from m4.utils.saveIFFInfo import SaveAdditionalInfo
from m4.utils import logger
import numpy as np
import pyfits
import os
from m4.utils.interferometer_converter import InterferometerConverter

class AnalyzerIFF(object):
    
    def __init__(self):
        self._indexingList= None
        self._cube= None
        self._reducedCube= None
        self._rec= None
        self._intMat= None

    def _ttData(self):
        split= os.path.split(self._h5Folder)
        self._tt= split[1]
        return self._tt
    
    @staticmethod
    def loadModalIFFInfoFromH5Folder(h5Folder):
        theObject= AnalyzerIFF()
        theObject._h5Folder= h5Folder
        a= SaveAdditionalInfo.loadIFFInfoModal(h5Folder)
        theObject._who= a[0]
        theObject._actsVector= a[1]
        theObject._cmdMatrix= a[2]
        theObject._indexingList= a[3]
        theObject._cmdAmplitude= a[4]
        theObject._nPushPull= a[5]
        return theObject
    
    @staticmethod
    def loadZonalIFFInfoFromH5Folder(h5Folder):
        theObject= AnalyzerIFF()
        theObject._h5Folder= h5Folder
        a= SaveAdditionalInfo.loadIFFInfoZonal(h5Folder)
        theObject._who= a[0]
        theObject._actsVector= a[1]
        theObject._cmdAmplitude= a[2]
        theObject._nPushPull= a[3]
        return theObject
        
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
        
        aa=np.arange(self._actsVector.shape[0])
        zipped= zip(aa, self._cmdAmplitude)
        for i, amp in zipped:
            for k in range(self._nPushPull):
                if self._indexingList is None:
                    mis = self._actsVector[i]
                    nomePos= 'pos%03d_pp%03d.h5' % (mis, k)               
                    nomeNeg= 'neg%03d_pp%03d.h5' % (mis, k)
                else:
                    p= self._nPushPull * i + k
                    where= self._indexReorganization()
                    n=where[p][0][0]
                    mis= k* self._indexingList.shape[1] + n
                
                    nomePos= 'pos%03d.h5' % mis                
                    nomeNeg= 'neg%03d.h5' % mis
                    
                imgPos= ic.from4D(os.path.join(self._h5Folder,nomePos))
                imgNeg= ic.from4D(os.path.join(self._h5Folder,nomeNeg))
                imgIF= (imgPos - imgNeg) / (2*amp)
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
                
        return 
    
    def getCube(self):
        return self._cube
    
    def getMasterMask(self):
        aa=np.sum(self._cube.mask.astype(int), axis=2)
        masterMask= np.zeros(aa.shape, dtype=np.bool)
        masterMask[np.where(aa>0)]= True
        return masterMask
    
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


  
    def ifRedux(self, masterROI, secondaryROI,
                solvePhaseAmbiguity= None, removeResidualTilt= None):
        return self._reducedCube
    
    def setAnalysisMaskFromMasterMask(self):
        self._analysisMask= self.getMasterMask()
    
    def _createInteractionMatrix(self):
        if self._analysisMask is None:
            self.setAnalysisMaskFromMasterMask()
        nActsInCube= self.getCube().shape[2] 
        pass
    
    def _createSurfaceReconstructor(self):
        pass
    
    def getInteractionMatrix(self):
        if self._intMat is None:
            self._createInteractionMatrix()
        return self._intMat

    
    def getReconstructor(self):
        if self._rec is None:
            self._createSurfaceReconstructor()
        return self._rec