'''
@author: cs
'''

import numpy as np
from m4.utils import saveIFFInfo
from m4.utils import logger
import copy
import os
from m4.type.commandHistory import CmdHistory

class IFFunctionsMaker(object):


    def __init__(self, device):
        self._device= device
        self._who= self._device._who
        self._nActs = self._device.nActs()

    
    
    def acq_IFFunctions(self, storeInFolder, modesVector, nPushPull, amplitude, 
                        cmdMatrix):
        
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
        '''
        indexingImput= copy.copy(modesVector)
        save= saveIFFInfo.SaveAdditionalInfo(storeInFolder)
        dove, tt= save._createFolderToStoreMeasurements()
        
        if np.count_nonzero(cmdMatrix - np.diag(np.diagonal(cmdMatrix))) == 0:
            print ('Misuro IF zonali')
            logger.log("Misura delle funzioni di influenza zonali", self._who, tt)
            save._saveIFFInfo(dove, self._who, amplitude,
                                indexingImput, nPushPull, cmdMatrix)
            
            
        else:
            print ('Misuro IF globali')
            logger.log("Misura delle funzioni di influenza modali", self._who, tt)
            
            cmdH= CmdHistory(self.device)
            
            
            
            
            
            
            
            
            