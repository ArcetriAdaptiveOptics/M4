'''
@author: cs
'''

import numpy as np
from m4.copies import C_saveIFFInfo
from m4.utils import logger
import copy
import os
from m4.type.commandHistory import CmdHistory

class IFFunctionsMaker(object):


    def __init__(self, device):
        self._device= device
        self._who= self._device._who
        self._nActs = self._device.nActs()
        
    
    def _pokeActuators(self, indexing, amplitude, 
                            pushOrPull, modeMatrix= None):
        self._device.pokeActs(indexing, amplitude, 
                             pushOrPull, modeMatrix=None)
        
    def _measureAndStoreH5(self, filename):
        #deve acquisire e mettere il file h5 nel posto giusto
        pass
    
    def _amplitudeReorganization(self, indexing, indexingList, 
                                 amplitude, nPushPull):
        where=[] 
        for i in indexing:           
            for j in range(nPushPull): 
                a=np.where(indexingList[j]==i)  
                where.append(a) 
         
         
        where=np.array(where)
        vect=np.zeros(amplitude.shape[0]*nPushPull) 
          
        for i in range(amplitude.shape[0]): 
            for k in range(nPushPull): 
                p= nPushPull * i + k
                indvect=where[p][0][0]+ indexing.shape[0] * k
                vect[indvect]=amplitude[i]
        
        return vect
  
    
    def acq_IFFunctions(self, storeInFolder, indexing, nPushPull, amplitude, 
                        cmdMatrix= None):
        
        '''
         arg:
             indexing= vettore dell'indice dei modi o degli attuatori
                        da applicare (numpy.array([]))
             nPushPull= numero di push pull consecutivi sull'attuatore
                         (int)
             Amplitude= vettore con l'ampiezza dei modi (numpy.array([]))
             cmdMatrix= matrice dei comandi dei modi 
                         (nActs x nModes)
        '''
        indexingImput= copy.copy(indexing)
        save= saveIFFInfo.SaveAdditionalInfo(storeInFolder)
        dove, tt= save._createFolderToStoreMeasurements()
        vecPushPull= np.array((1,-1)) 
        
        
        if cmdMatrix is None:
            print ('Misuro IF zonali')
            logger.log("Misura delle funzioni di influenza zonali", self._who, tt)
            save._saveIFFInfoZonal(dove, self._who, amplitude,
                                indexingImput, nPushPull)
            
            zipped=zip(indexing, amplitude)
            for ind, amp in zipped:
                for j in range(nPushPull):
                    self._pokeActuators(ind, amp,
                                         vecPushPull[0])
                    nome='pos%03d_pp%03d' % (ind, j)
                    self._measureAndStoreH5(os.path.join(dove, nome))
                    print('misura positiva attuatore %03d, Push %03d' %(ind, j))
                    self._pokeActuators(ind, amp, 
                                              vecPushPull[1])
                    nome='neg%03d_pp%03d' % (ind, j)
                    self._measureAndStoreH5(os.path.join(dove, nome))
                    print('misura negativa attuatore %03d, Push %03d' %(ind, j))
                        
        
        else:
            print ('Misuro IF globali')
            logger.log("Misura delle funzioni di influenza modali", self._who, tt)
            
            cmdH= CmdHistory(self.device)
            matrixToApply, indexingList= cmdH.shuffleCmdHistoryMaker(
                                            indexing, nPushPull, cmdMatrix)
            cmdH.saveCmdHistory()
             
            randomIndexing= np.array(indexingList)     
            save._saveIFFInfoModal(dove, self._who, amplitude,
                                indexingImput, nPushPull, 
                                randomIndexing, cmdMatrix) 
              
            aa= np.arange(matrixToApply.shape[1])
            reorganizedAmplitude= self._amplitudeReorganization(indexingImput, 
                                                indexingList, amplitude, nPushPull)
            zipped= zip(aa, reorganizedAmplitude)
            
            for ind, amp in zipped:
                
                self._pokeActuators(ind, amp,
                                         vecPushPull[0], matrixToApply)
                nome='pos%03d' % ind
                self._measureAndStoreH5(os.path.join(dove, nome))
                print('misura positiva %03d' %ind)
                self._pokeActuators(ind, amp, 
                                        vecPushPull[1], matrixToApply)
                nome='neg%03d' % ind
                self._measureAndStoreH5(os.path.join(dove, nome))
                print('misura negativa %03d' %ind)
                
            
        return dove

    
            
