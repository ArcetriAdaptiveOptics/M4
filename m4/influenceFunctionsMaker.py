'''
@author: cs
'''

import numpy as np
from m4.ground import trackingNumberFolder
from m4.ground.configuration import Configuration
from m4.ground import logger
import pyfits
import copy
import os
import h5py
from m4.type.commandHistory import CmdHistory
from m4.ground import objectFromFitsFileName


class IFFunctionsMaker(object):

    def __init__(self, device):
        self._device= device
        self._who= self._device._who
        self._nActs = self._device.nActs()
        
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "IFFunctions")
    
    
    def acq_IFFunctions(self, modesVectorFitsFileName, nPushPull, amplitudeFitsFileName, 
                        cmdMatrixFitsFileName, shuffle=None):
        
        amplitude, modesVector, cmdMatrix= objectFromFitsFileName.readObjectFitsFileName(amplitudeFitsFileName, modesVectorFitsFileName,
                                                                        cmdMatrixFitsFileName)
        self._nPushPull= nPushPull
        
        self._modesVectorTag= modesVectorFitsFileName
        self._amplitudeTag= amplitudeFitsFileName
        self._cmdMatrixTag= cmdMatrixFitsFileName
        '''
         arg:
             modesVectorTag= fits file name (modes.fits)
             nPushPull= numero di push pull consecutivi sull'attuatore
                         (int)
             amplitudeTag= fits file name (amp.fits)
             cmdMatrixTag= fits file name (modalBase.fits)
             shuffle= se non indicato viene creata la matrice della storia dei comandi ordinata
             
        return:
                tracking number delle misure effettuate
        '''
        storeInFolder= self._storageFolder()
        indexingInput= copy.copy(modesVector)
        save= trackingNumberFolder.TtFolder(storeInFolder)
        dove, tt= save._createFolderToStoreMeasurements()
        
        diagonalMat= self._diagonalControll(cmdMatrix)
        if np.count_nonzero(diagonalMat - np.diag(np.diagonal(diagonalMat))) == 0:
            print ('Measure of zonal IF')
            logger.log("Measurement of zonal influence functions", self._who, tt)
            
        else:
            print ('Measure of global IF')
            logger.log("Measurement of modal influence functions", self._who, tt)
        
                
        cmdH= CmdHistory(self._device)
        
        if shuffle is None:
            commandHistoryMatrixToApply, self._tt_cmdH= cmdH.tidyCommandHistoryMaker(modesVector,
                                                                  amplitude,
                                                                  cmdMatrix, 
                                                                  nPushPull)
        else:
            commandHistoryMatrixToApply, self._tt_cmdH= cmdH.shuffleCommandHistoryMaker(modesVector,
                                                                  amplitude,
                                                                  cmdMatrix, 
                                                                  nPushPull)
        self._indexingList= cmdH.getIndexingList()
        self._saveInfoAsFits(dove)
         
         
        self._applyToDM() 
        
        return tt
           
           
    def _diagonalControll(self, matrix):
        v= np.zeros((self._nActs,1))
        reps= matrix.shape[0] - matrix.shape[1]
        vects= np.tile(v, reps)
        newMatrix= np.hstack((matrix, vects))
        
        return newMatrix
        
     
    def _applyToDM(self):
        pass
    
    def _testIFFunctions_createCube25fromFileFitsMeasure(self):
        fold='/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/OIM_25modes.fits' 
        hduList= pyfits.open(fold) 
        cube50= hduList[0].data
        
        imaList=[]
        maskList=[]
        for i in range(cube50.shape[0]):
            if i%2 == 0:
                imaList.append(cube50[i])
            else:
                maskList.append(cube50[i])
                
        cube25= None
        zipped= zip(imaList, maskList)      
        for ima, mask in zipped:
            immagine=np.ma.masked_array(ima, mask= mask)
            if cube25 is None:
                cube25= immagine
            else:
                cube25= np.ma.dstack((cube25, immagine))
        return cube25
     
            
    def _saveInfoAsFits(self, folder):
        fitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.Header()
        header['NPUSHPUL']= self._nPushPull
        header['WHO']= self._who
        header['TT_CMDH']= self._tt_cmdH
        header['MODEVECT']= self._modesVectorTag
        header['CMDMAT']= self._cmdMatrixTag
        header['AMP']= self._amplitudeTag
        pyfits.writeto(fitsFileName, self._indexingList, header)
        
    def _saveInfoAsH5(self, folder):
        fitsFileName= os.path.join(folder, 'info.h5')
        hf = h5py.File(fitsFileName, 'w')
        hf.create_dataset('dataset_1', data=self._indexingList)
        hf.attrs['MODEVECT']= self._modesVectorTag
        hf.attrs['CMDMAT']= self._cmdMatrixTag
        hf.attrs['AMP']= self._amplitudeTag
        hf.attrs['NPUSHPUL']= self._nPushPull
        hf.attrs['WHO']= self._who
        hf.attrs['TT_CMDH']= self._tt_cmdH
        hf.close()

            
    @staticmethod
    def loadInfoFromFits(folder):
        additionalInfoFitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.getheader(additionalInfoFitsFileName)
        hduList= pyfits.open(additionalInfoFitsFileName)
        actsVectorTag= header['MODEVECT']
        cmdMatrixTag= header['CMDMAT']
        cmdAmplTag= header['AMP']
        indexingList= hduList[0].data
        who= header['WHO']
        tt_cmdH= header['TT_CMDH']
        try:
            nPushPull= header['NPUSHPUL']
        except KeyError:
            nPushPull= 1
            
        cmdAmpl, actsVector, cmdMatrix= objectFromFitsFileName.readObjectFitsFileName(cmdAmplTag, actsVectorTag, cmdMatrixTag)
        return (who, tt_cmdH, actsVector, cmdMatrix, cmdAmpl, nPushPull, indexingList)
    
    @staticmethod
    def loadInfoFromH5(tt):
        storeInFolder= IFFunctionsMaker._storageFolder()
        folder= os.path.join(storeInFolder, tt)
        fileName= os.path.join(folder, 'info.h5')
        hf = h5py.File(fileName, 'r')
        hf.keys()
        data1= hf.get('dataset_1')
        actsVectorTag= hf.attrs['MODEVECT']
        cmdMatrixTag= hf.attrs['CMDMAT']
        cmdAmplTag= hf.attrs['AMP']
        indexingList= np.array(data1)
        nPushPull= hf.attrs['NPUSHPUL']
        who= hf.attrs['WHO']
        tt_cmdH= hf.attrs['TT_CMDH']
        
        cmdAmpl, actsVector, cmdMatrix= objectFromFitsFileName.readObjectFitsFileName(cmdAmplTag, actsVectorTag, cmdMatrixTag)
        return (who, tt_cmdH, actsVector, cmdMatrix, cmdAmpl, nPushPull, indexingList)
     
            
            
            
            
            