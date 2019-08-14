'''
@author: cs
'''
import m4 as _m4
import numpy as np

class IFFunctions(object):


    def __init__(self, device):
        self._device= device
        self._nActs = self._device.nActs()
        
    
    def _pokeActuators(self, indexing, amplitude, 
                            pushOrPull, modeMatrix= None):
        self._device.pokeActs(indexing, amplitude, 
                             pushOrPull, modeMatrix=None)
  
    
    def pushPull(self, indexing, nPushPull, amplitude, 
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
        vecPushPull= np.array((1,-1)) 
        
        if cmdMatrix is None:
            print ('Misuro IF zonali')
            
            for ind in indexing:
                for j in range(nPushPull):
                    self._pokeActuators(ind,amplitude,
                                         vecPushPull[0])
                    self._pokeActuators(ind, amplitude, 
                                              vecPushPull[1])
                  
        
        else:
            print ('Misuro IF globali')
            
            nFrame= indexing.size * nPushPull
            matrixToApply= np.zeros((self._nActs,nFrame))
            
            cmdList=[]
            for i in indexing:
                cmd= cmdMatrix[:,i]
                cmdList.append(cmd)
                
            for j in range(nPushPull):
                     
                for i in range(len(cmdList)):
                    k= len(cmdList)*j + i
                    matrixToApply.T[k]=cmdList[i]
                    
                  
            for ind in indexing:
                for j in range(nPushPull):
                    self._pokeActuators(ind, amplitude,
                                         vecPushPull[0], matrixToApply)
                    self._pokeActuators(ind, amplitude, 
                                        vecPushPull[1], matrixToApply)
                
            
            return matrixToApply
        
                
    def acqIntMat(self):
        pass
    
    
    def ifRedux(self):
        pass
    
    
    def buildIntMat(self):
        pass
    
    
    def buildRecMat(self):
        pass
    
    
            
