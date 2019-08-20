'''
@author: cs
'''
from m4.utils.saveInfo import SaveAdditionalInfo
from m4.utils import logging
import numpy as np
import os
from m4.utils.interferometer_converter import InterferometerConverter

class AnalyzerIFF(object):
    
    def __init__(self, device):
        self._device= device[0]
        self._who= device[1]
        self._nActs = self._device.nActs()
        self._cube= None
        self._rec= None
        self._intMat= None
        
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
        logging.log('Creo il cubo delle IFF per', self._who)
        
        for i in range(self._actsVector.shape[0]):
            for k in range(self._nPushPull):
                p= self._nPushPull * i + k
                where= self._indexReorganization()
                n=where[p][0][0]
                
                nomePos= 'pos%03d.h5' % n                
                nomeNeg= 'neg%03d.h5' % n
                imgPos= ic.from4D(os.path.join(self._h5Folder,nomePos))
                imgNeg= ic.from4D(os.path.join(self._h5Folder,nomeNeg))
                imgIF= (imgPos - imgNeg) / (2*self._cmdAmplitude)
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
    
    @staticmethod
    def fromH5Folder(device, h5Folder):
        theObject= AnalyzerIFF(device)
        theObject._h5Folder= h5Folder
        a= SaveAdditionalInfo.loadAdditionalInfo(h5Folder)
        theObject._who= a[0]
        theObject._actsVector= a[1]
        theObject._cmdMatrix= a[2]
        theObject._indexingList= a[3]
        theObject._cmdAmplitude= a[4]
        theObject._nPushPull= a[5]
        return theObject
    
    def _createInteractionMatrix(self):
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