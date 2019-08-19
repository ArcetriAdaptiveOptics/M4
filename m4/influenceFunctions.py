'''
@author: cs
'''

import numpy as np
from m4.utils import saveInfo
from m4.utils import logging
import copy

class IFFunctions(object):


    def __init__(self, device):
        self._device= device
        self._nActs = self._device.nActs()
        
    
    def _pokeActuators(self, indexing, amplitude, 
                            pushOrPull, modeMatrix= None):
        self._device.pokeActs(indexing, amplitude, 
                             pushOrPull, modeMatrix=None)
  
    
    def pushPull(self, storeInFolder, indexing, nPushPull, amplitude, 
                 cmdMatrix= None):
        
        '''
         arg:
             indexing= vettore dell'indice dei modi o degli attuatori
                        da applicare (numpy.array([]))
             nPushPull= numero di push pull consecutivi sull'attuatore
                         (int)
             Amplitude= ampiezza del modo
             cmdMatrix= matrice dei comandi dei modi 
                         (nActs x nModes)
        '''
        indexingImput= copy.copy(indexing)
        save= saveInfo.SaveAdditionalInfo(storeInFolder)
        vecPushPull= np.array((1,-1)) 
        
        
        if cmdMatrix is None:
            print ('Misuro IF zonali')
            logging.log("Misuro le funzioni di influenza", "zonali")
            save._saveIFFInfo(storeInFolder, amplitude,
                                indexingImput, nPushPull)
            
            for ind in indexing:
                for j in range(nPushPull):
                    self._pokeActuators(ind,amplitude,
                                         vecPushPull[0])
                    self._pokeActuators(ind, amplitude, 
                                              vecPushPull[1])
                  
        
        else:
            print ('Misuro IF globali')
            logging.log("Misuro le funzioni di influenza", "globali")
            
            nFrame= indexing.size * nPushPull
            matrixToApply= np.zeros((self._nActs,nFrame))
             
            indexingList=[]     
            for j in range(nPushPull):
                np.random.shuffle(indexing)
                indexingList.append(list(indexing))
                 
                cmdList=[]
                for i in indexing:
                    cmd= cmdMatrix[:,i]
                    cmdList.append(cmd)
                
                for i in range(len(cmdList)):
                    k= len(cmdList)*j + i
                    matrixToApply.T[k]=cmdList[i]
             
            randomIndexing= np.array(indexingList)     
            save._saveIFFInfo(storeInFolder, amplitude,
                                indexingImput, nPushPull, 
                                randomIndexing, cmdMatrix)   
                  
            for ind in indexing:
                for j in range(nPushPull):
                    self._pokeActuators(ind, amplitude,
                                         vecPushPull[0], matrixToApply)
                    self._pokeActuators(ind, amplitude, 
                                        vecPushPull[1], matrixToApply)
                
            
            return matrixToApply, indexingList
        
                
    def acqIntMat(self):
        pass
    
    
    def ifRedux(self):
        pass
    
    
    def buildIntMat(self):
        pass
    
    
    def buildRecMat(self):
        pass
    
    
            
