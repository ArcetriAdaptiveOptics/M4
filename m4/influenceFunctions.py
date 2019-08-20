'''
@author: cs
'''

import numpy as np
from m4.utils import saveInfo
from m4.utils import logging
import copy
import os

class IFFunctions(object):


    def __init__(self, device):
        self._device= device[0]
        self._who= device[1]
        self._nActs = self._device.nActs()
        
    
    def _pokeActuators(self, indexing, amplitude, 
                            pushOrPull, modeMatrix= None):
        self._device.pokeActs(indexing, amplitude, 
                             pushOrPull, modeMatrix=None)
        
    def _measureAndStoreH5(self, filename):
        #deve acquisire e mettere il file h5 nel posto giusto
        pass
  
    
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
        dove= save._createFolderToStoreMeasurements()
        vecPushPull= np.array((1,-1)) 
        
        
        if cmdMatrix is None:
            print ('Misuro IF zonali')
            logging.log("Misuro le funzioni di influenza", "zonali")
            save._saveIFFInfo(dove, self._who, amplitude,
                                indexingImput, nPushPull)
            
            for ind in indexing:
                for j in range(nPushPull):
                    self._pokeActuators(ind,amplitude,
                                         vecPushPull[0])
                    nome='pos%03d_pp%03d' % (ind, j)
                    self._measureAndStoreH5(os.path.join(dove, nome))
                    print('misura positiva %03d, Push %03d' %(ind, j))
                    self._pokeActuators(ind, amplitude, 
                                              vecPushPull[1])
                    nome='neg%03d_pp%03d' % (ind, j)
                    self._measureAndStoreH5(os.path.join(dove, nome))
                    print('misura negativa %03d, Push %03d' %(ind, j))
                  
        
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
            save._saveIFFInfo(dove, self._who, amplitude,
                                indexingImput, nPushPull, 
                                randomIndexing, cmdMatrix)   
                  
            for ind in range(matrixToApply.shape[1]):
                
                self._pokeActuators(ind, amplitude,
                                         vecPushPull[0], matrixToApply)
                nome='pos%03d' % i
                self._measureAndStoreH5(os.path.join(dove, nome))
                print('misura positiva %03d' %ind)
                self._pokeActuators(ind, amplitude, 
                                        vecPushPull[1], matrixToApply)
                nome='neg%03d' % i
                self._measureAndStoreH5(os.path.join(dove, nome))
                print('misura negativa %03d' %ind)
                
            
            return matrixToApply, indexingList

    
    
    def ifRedux(self):
        pass

    
            
