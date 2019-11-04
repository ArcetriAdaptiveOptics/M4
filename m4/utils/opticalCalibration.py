'''
@author: cselmi
'''

from m4.ground.configuration import Configuration
from m4.utils.zernikeOnM4 import ZernikeOnM4
from m4.ground import trackingNumberFolder
from m4.ground import logger
import numpy as np
import pyfits
import os


class Opt_Calibration():
    
    def __init__(self):
        self._zOnM4= ZernikeOnM4()
        self._rec= None
        self._intMat= None
        
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "Calibration")
        
        
    def measureCalibrationMatrix(self, who, commandAmpVector, nPushPull):
        self._nPushPull= nPushPull
        self._commandAmpVector= commandAmpVector
        '''
            arg:
                who= numero che indica l'elemento ottico su cui svolgere la calibrazione
                commandAmpVector= vettore contenente l'ampiezza dei comandi da dare ai gradi di libertà da calibrare
        '''
        storeInFolder= self._storageFolder()
        save= trackingNumberFolder.TtFolder(storeInFolder)
        dove, self._tt= save._createFolderToStoreMeasurements()
        logger.log('Misure di', 'calibrazione', 'con tt=', self._tt)
        
        self._commandMatrix= self._createCommandMatrix(who, self._commandAmpVector, self._nPushPull)
        self._saveCommandMatrixAsFits(dove)
        
        self._applyCommandMatrix(self._commandMatrix)
        return self._tt
    
    def analyzerCalibrationMeasurement(self, tt, maskIndex):
        '''
        arg:
             maskIndex= int dell'indice della maschera relativo allo specchio di riferimento
        '''
        a= Opt_Calibration.loadCommandMatrixFromFits(tt)
        a.createCube(tt)
        cube= a.getCube()
        ima= cube[:,:,0]
        from m4.utils.roi import ROI 
        r= ROI()
        roi= r.ROIonAlignmentImage(ima)  
        mask= roi[maskIndex]
        dove= os.path.join(self._storageFolder(), tt)
        self._saveMask(dove, mask)
        self._intMat= a.getInteractionMatrix(mask)
        self._rec= a.getReconstructor(mask)
        self._saveIntMatAndRec(dove)
        return self._intMat, self._rec
        

    def _applyCommandMatrix(self, commandMatrix):
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
            self._rows= Configuration.ParabolDoF + Configuration.RMDof
        elif who==1:
            self._who= 'PAR'
            self._rows= Configuration.ParabolDoF
        elif who==2:
            self._who= 'RM'
            self._rows= Configuration.RMDof
        elif who==3:
            self._who= 'M4'
            self._rows= Configuration.M4Dof
        else:
            raise OSError('Who= %s doesnt exists' % who)  
        
        logger.log('Creazione delle matrice', 'dei comandi', 'per', self._who) 
        self._commandMatrix= self._createCommandHistoryMatrix(self._rows, commandAmpVector, nPushPull)   
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
      
    
    def _saveCommandMatrixAsFits(self, dove):   
        fitsFileName= os.path.join(dove, 'CommandMatrix.fits')
        header= pyfits.Header()
        header['NPUSHPUL']= self._nPushPull
        header['WHO']= self._who
        pyfits.writeto(fitsFileName, self._commandAmpVector, header)
        pyfits.append(fitsFileName, self._commandMatrix, header)
        
    @staticmethod
    def loadCommandMatrixFromFits(tt):
        theObject= Opt_Calibration()
        theObject._tt= tt
        dove= os.path.join(theObject._storageFolder(), tt)
        file= os.path.join(dove, 'CommandMatrix.fits')
        header= pyfits.getheader(file)
        hduList= pyfits.open(file)
        theObject._who= header['WHO']
        theObject._nPushPull= header['NPUSHPUL']
        theObject._commandAmpVector= hduList[0].data
        theObject._commandMatrix= hduList[1].data
        return theObject
    
    def _testCalibration_createCubeMeasurefromFileFitsMeasure(self):
        cubeMeasure= None 
        fold='/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/MixingIntMat/20190930_162714' 
        for i in range(5):
            name= 'Frame_%04d.fits' %i
            file= os.path.join(fold, name)
            hduList= pyfits.open(file)
            aa= hduList[0].data
            mode= np.ma.masked_array(aa[0,:,:], mask= np.invert(aa[1,:,:].astype(bool)))
            if cubeMeasure is None:
                cubeMeasure= mode
            else:
                cubeMeasure= np.ma.dstack((cubeMeasure, mode))
        return cubeMeasure
     
    def createCube(self, tt):
        logger.log('Creazione del', 'cubo', 'relativo a', self._tt)
        cubeFromMeasure= self._testCalibration_createCubeMeasurefromFileFitsMeasure() 
        for i in range(cubeFromMeasure.shape[2]):
            cubeFromMeasure[:,:,i]= cubeFromMeasure[:,:,i] / self._commandAmpVector[i]
        self._cube= cubeFromMeasure
        return 
    
    def getCube(self):
        return self._cube  
      
    def _createInteractionMatrixAndReconstructor(self, mask):
        coefList=[]
        for i in range(self._cube.shape[2]):
            ima= np.ma.masked_array(self._cube[:,:,i], mask=mask)
            coef, mat= self._zOnM4.zernikeFit(ima, np.arange(2,11))
            #z= np.array([2,3,4,7,8])
            z= np.array([0,1,2,5,6])
            finalCoef= np.zeros(z.shape[0])
            aa= np.arange(finalCoef.shape[0])
            zipped=zip(aa,z)
            for i,j in zipped:
                finalCoef[i]=coef[j]
            self._mat=mat
            coefList.append(finalCoef)
        
        self._intMat= np.zeros((coefList[0].shape[0], self._cube.shape[2]))
        for j in range(self._cube.shape[2]):
            self._intMat.T[j]= coefList[j]
            
        self._rec= np.linalg.pinv(self._intMat)
        
    
    def getInteractionMatrix(self, mask):
        if self._intMat is None:
            self._createInteractionMatrixAndReconstructor(mask)
        return self._intMat
    
    def getReconstructor(self, mask):
        if self._rec is None:
            self._createInteractionMatrixAndReconstructor(mask)
        return self._rec
    
    def _saveIntMatAndRec(self, dove):
        fitsFileName= os.path.join(dove, 'InteractionMatrix.fits')
        pyfits.writeto(fitsFileName, self._intMat)
        fitsFileName= os.path.join(dove, 'Reconstructor.fits')
        pyfits.writeto(fitsFileName, self._rec)
        
    def _saveMask(self, dove, mask):
        fitsFileName= os.path.join(dove, 'Mask.fits')
        pyfits.writeto(fitsFileName, mask.astype(int))
        
        
        
        
        
        
        