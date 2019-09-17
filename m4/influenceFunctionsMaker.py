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
from m4.type.commandHistory import CmdHistory

class IFFunctionsMaker(object):


    def __init__(self, device):
        self._device= device
        self._who= self._device._who
        self._nActs = self._device.nActs()
        
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "IFFunctions")

    
    
    def acq_IFFunctions(self, modesVector, nPushPull, amplitude, 
                        cmdMatrix, shuffle=None):
        
        '''
         arg:
             modesVector= vettore dell'indice dei modi o degli attuatori
                        da applicare (numpy.array([]))
             nPushPull= numero di push pull consecutivi sull'attuatore
                         (int)
             Amplitude= vettore con l'ampiezza dei modi (numpy.array([]))
             cmdMatrix= matrice dei comandi dei modi 
                         (nActs x nModes)
                         matrice diagonale nel caso di comandi zonali
             shuffle= se non indicato viene creata la matrice della storia dei comandi ordinata
             
        return:
                tracking number delle misure effettuate
        '''
        storeInFolder= self._storageFolder()
        indexingImput= copy.copy(modesVector)
        save= trackingNumberFolder.TtFolder(storeInFolder)
        dove, tt= save._createFolderToStoreMeasurements()
        
        diagonalMat= self._diagonalControll(cmdMatrix)
        if np.count_nonzero(diagonalMat - np.diag(np.diagonal(diagonalMat))) == 0:
            print ('Misuro IF zonali')
            logger.log("Misura delle funzioni di influenza zonali", self._who, tt)
            
        else:
            print ('Misuro IF globali')
            logger.log("Misura delle funzioni di influenza modali", self._who, tt)
        
                
        cmdH= CmdHistory(self._device)
        
        if shuffle is None:
            commandHistoryMatrixToApply, tt_cmdH= cmdH.tidyCommandHistoryMaker(modesVector,
                                                                  amplitude,
                                                                  cmdMatrix, 
                                                                  nPushPull)
        else:
            commandHistoryMatrixToApply, tt_cmdH= cmdH.shuffleCommandHistoryMaker(modesVector,
                                                                  amplitude,
                                                                  cmdMatrix, 
                                                                  nPushPull)
        self._indexingList= cmdH.getIndexingList()
        self._saveInfo(dove, self._who, tt_cmdH, amplitude, indexingImput,
                        nPushPull, cmdMatrix, self._indexingList)
         
         
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
     
            
    def _saveInfo(self, folder, who, tt_cmdH, amplitude,vectorOfActuatorsOrModes,
                     nPushPull, cmdMatrix, indexingList):
            
        fitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.Header()
        header['NPUSHPUL']= nPushPull
        header['WHO']= who
        header['TT_CMDH']= tt_cmdH
        pyfits.writeto(fitsFileName, vectorOfActuatorsOrModes, header)
        pyfits.append(fitsFileName, cmdMatrix, header)
        pyfits.append(fitsFileName, amplitude, header)
        pyfits.append(fitsFileName, indexingList, header)

            
    @staticmethod
    def loadInfo(folder):
        additionalInfoFitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.getheader(additionalInfoFitsFileName)
        hduList= pyfits.open(additionalInfoFitsFileName)
        actsVector= hduList[0].data
        cmdMatrix= hduList[1].data
        cmdAmpl= hduList[2].data
        indexingList= hduList[3].data
        who= header['WHO']
        tt_cmdH= header['TT_CMDH']
        try:
            nPushPull= header['NPUSHPUL']
        except KeyError:
            nPushPull= 1
            
        return (who, tt_cmdH, actsVector, cmdMatrix, cmdAmpl, nPushPull, indexingList)
     
            
            
            
            
            