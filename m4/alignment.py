'''
@author: cselmi
'''

from m4.ground.configuration import Configuration
from m4.ground.zernikeGenerator import ZernikeGenerator
from m4.ground import trackingNumberFolder
from m4.ground import logger
import numpy as np
import pyfits
import os


class Alignment():
    
    def __init__(self):
        self._pupilXYRadius= Configuration.ParabolaPupilXYRadius
        self._zg= ZernikeGenerator(2*self._pupilXYRadius[2])
        self._rec= None
        self._intMat= None
        
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "Alignment")
        
        
    def measureCalibrationMatrix(self, who, commandAmpVector, nPushPull):
        self._nPushPull= nPushPull
        self._commandAmpVector= commandAmpVector
        '''
            arg:
                mask= measurement mask (segment mask or RM mask)
                commandAmpVector= vettore contenente l'ampiezza dei comandi da dare ai gradi di libertà da calibrare
        '''
        storeInFolder= self._storageFolder()
        save= trackingNumberFolder.TtFolder(storeInFolder)
        dove, self._tt= save._createFolderToStoreMeasurements()
        
        self._commandMatrix= self._createCommandMatrix(who, self._commandAmpVector, self._nPushPull)
        self.saveCommandMatrixAsFits(dove)
        
        self._applyCommandMatrixToDM(self._commandMatrix)
        return self._tt
        

    def _applyCommandMatrixToDM(self, commandMatrix):
        #deve applicare la matrice e salvare gli interferogrammi
        pass   
    
    
    def _createCommandMatrix(self, who, commandAmpVector, nPushPull):
        '''
            arg:
                who=
                    0 per mixing
                    1 per parabola
                    2 per specchio di riferimento
                    3 per specchio deformabile
        '''
        if who==0:
            self._who= 'PAR + RM'
            rows= Configuration.ParabolDoF + Configuration.RMDof
            self._commandMatrix= self._createCommandHistoryMatrix(rows, commandAmpVector, nPushPull)
        elif who==1:
            self._who= 'PAR'
            rows= Configuration.ParabolDoF
            self._commandMatrix= self._createCommandHistoryMatrix(rows, commandAmpVector, nPushPull)
        elif who==2:
            self._who= 'RM'
            rows= Configuration.RMDof
            self._commandMatrix= self._createCommandHistoryMatrix(rows, commandAmpVector, nPushPull)
        elif who==3:
            self._who= 'M4'
            rows= Configuration.M4Dof
            self._commandMatrix= self._createCommandHistoryMatrix(rows, commandAmpVector, nPushPull)
        else:
            raise OSError('Who= %s doesnt exists' % who)  
        
        logger.log('Creazione delle matrice', 'dei comandi', 'per', self._who)    
        return self._commandMatrix
               
    
    def _createCommandHistoryMatrix(self, rows, commandAmpVector, nPushPull):
        '''
            crea la matrice dei comandi usando come righe i gradi di libertà scelti in row
        '''
        vecPushPull= np.array((1,-1)) 
        rows= rows
        columns= commandAmpVector.shape[0]*nPushPull*vecPushPull.shape[0]
        commandMatrix= np.zeros((rows, columns))
        
        for k in range(nPushPull):
            for i in range(commandAmpVector.shape[0]):
                j= (commandAmpVector.shape[0]*vecPushPull.shape[0])*k +2*i
                commandMatrix[i,j]= commandAmpVector[i]*vecPushPull[0]
                commandMatrix[i,j+1]= commandAmpVector[i]*vecPushPull[1]   
        return commandMatrix 
      
    
    def saveCommandMatrixAsFits(self, dove):   
        fitsFileName= os.path.join(dove, 'CommandMatrix.fits')
        header= pyfits.Header()
        header['NPUSHPUL']= self._nPushPull
        header['WHO']= self._who
        pyfits.writeto(fitsFileName, self._commandAmpVector, header)
        pyfits.append(fitsFileName, self._commandMatrix, header)
        
    @staticmethod
    def loadCommandMatrixFromFits(tt):
        theObject= Alignment()
        theObject._tt= tt
        dove= os.path.join(theObject._storageFolder(), tt)
        file= os.path.join(dove, 'CommandMatrix.fits')
        header= pyfits.getheader(file)
        hduList= pyfits.open(file)
        theObject._who= header['WHO']
        theObject._nPushPull= header['NPUSHPUL']
        theObject._commandAmpVector= hduList[0]
        theObject._commandMatrix= hduList[1]
        return theObject
    
    def _testAlignment_createCubeMeasurefromFileFitsMeasure(self):
        pass
     
    def createCube(self):
        totalCube=None 
        
        self._cube= totalCube               
        return self._cube  
      
    def _createInteractionMatrixAndReconstructor(self):
        
        pass
    
    def getInteractionMatrix(self):
        if self._intMat is None:
            self._createInteractionMatrixAndReconstructor()()
        return self._intMat
    
    def getReconstructor(self):
        if self._rec is None:
            self._createInteractionMatrixAndReconstructor()()
        return self._rec
       
        
        
        
        
        
        
        
        