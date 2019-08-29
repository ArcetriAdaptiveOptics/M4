'''
@author: cs
'''

from m4.utils import trackingNumberFolder
from m4.utils.configuration import Configuration
import numpy as np
from m4.utils import logger
import os 
import copy
import pyfits


class CmdHistory(object):
    
    def __init__(self, device):
        self._device= device
        self._who= self._device._who 
        self._nActs = self._device.nActs()
        self._modeVector= None 
        self._nPushPull= None
        self._cmdMatrix= None 
        self._cmdHToApply= None
        self._ampVect= None 
        
    def getCommandHistory(self):
        return self._cmdHToApply
        
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "CommandHistory")
        

    '''
         arg:
             modesVector= vettore dell'indice dei modi o degli attuatori
                        da applicare (numpy.array([]))
             nPushPull= numero di push pull consecutivi sull'attuatore
                         (int)
             ampVector= vettore con l'ampiezza dei modi (numpy.array([]))
             cmdMatrix= matrice dei comandi dei modi 
                         (nActs x nModes)
                         matrice diagonale nel caso di comandi zonali
    '''
        
    def tidyCommandHistoryMaker(self, modeVector, ampVector,
                                        cmdMatrix, nPushPull):
        self._ampVect= ampVector
        cmdHistory= self._tidyCmdHistory(modeVector, nPushPull, cmdMatrix)
        aa= np.arange(self._cmdHistory.shape[1])
        bb= np.tile(ampVector, nPushPull)
        zipped= zip(aa, bb)
        matrixToApply= self._cmdHistoryToApply(zipped)
        self._cmdHToApply= matrixToApply
        tt= self.save()
        logger.log('Creazione della commandHistoryMatrix', 'ordinata', tt)
        print(tt)
        
        return matrixToApply
    
    
    def shuffleCommandHistoryMaker(self, modeVector, ampVector, 
                                        cmdMatrix, nPushPull):
        self._ampVect= ampVector
        cmdHistory, indexingList= self._shuffleCmdHistory(
                                        modeVector, nPushPull, cmdMatrix)
        zipped= self._zippedAmplitude(ampVector)
        matrixToApply= self._cmdHistoryToApply(zipped)
        self._cmdHToApply=matrixToApply
        tt= self.save()
        logger.log('Creazione della commandHistoryMatrix', 'casuale', tt)
        print(tt)
        
        return matrixToApply
        
    
    def _shuffleCmdHistory(self, modeVector, nPushPull, cmdMatrix):
        self._modeVector= copy.copy(modeVector)
        self._nPushPull= nPushPull
        self._cmdMatrix= cmdMatrix
        
        nFrame= modeVector.size * nPushPull
        matrixToApply= np.zeros((self._nActs,nFrame))
             
        indexingList=[]     
        for j in range(nPushPull):
            np.random.shuffle(modeVector)
            indexingList.append(list(modeVector))
                 
            cmdList=[]
            for i in modeVector:
                cmd= cmdMatrix[:,i]
                cmdList.append(cmd)
                
            for i in range(len(cmdList)):
                k= len(cmdList)*j + i
                matrixToApply.T[k]=cmdList[i]
                
        self._cmdHistory= matrixToApply
        self._indexingList= np.array(indexingList)
                
        return self._cmdHistory, self._indexingList
    
    
    def _tidyCmdHistory(self, modeVector, nPushPull, cmdMatrix):
        self._modeVector= copy.copy(modeVector)
        self._nPushPull= nPushPull
        self._cmdMatrix= cmdMatrix
        self._indexingList= np.tile(modeVector, nPushPull)
        
        nFrame= modeVector.size * nPushPull
        matrixToApply= np.zeros((self._nActs,nFrame))
        
        cmdList=[]
        for i in modeVector:
            cmd= cmdMatrix[:,i]
            cmdList.append(cmd)
        
        for j in range(nPushPull):
            for i in range(len(cmdList)):
                k= len(cmdList)*j + i
                matrixToApply.T[k]=cmdList[i]
        
        self._cmdHistory= matrixToApply
        
        return self._cmdHistory
    
               
    def _amplitudeReorganization(self, indexingInput, indexingList, 
                                 amplitude, nPushPull):
        where=[] 
        for i in indexingInput:           
            for j in range(nPushPull): 
                a=np.where(indexingList[j]==i)  
                where.append(a) 
         
         
        where=np.array(where)
        vect=np.zeros(amplitude.shape[0]*nPushPull) 
          
        for i in range(amplitude.shape[0]): 
            for k in range(nPushPull): 
                p= nPushPull * i + k
                indvect=where[p][0][0]+ indexingInput.shape[0] * k
                vect[indvect]=amplitude[i]
        
        return vect
    
    
    def _zippedAmplitude(self, amplitude):
        aa= np.arange(self._cmdHistory.shape[1])
        reorganizedAmplitude= self._amplitudeReorganization(self._modeVector, 
                                                self._indexingList, 
                                                amplitude, self._nPushPull)
        #zipped= np.dstack((aa, reorganizedAmplitude))
        zipped= zip(aa, reorganizedAmplitude)
        return zipped
    
    def _cmdHistoryToApply(self, zipped):
        matrixWithAmp= self._cmdHistory
        for i, amp in zipped:
            matrixWithAmp.T[i]= matrixWithAmp.T[i]* amp
            
        vecPushPull= np.array((1,-1)) 
        matrixToApply= np.zeros((self._nActs,
                            self._cmdHistory.shape[1]*vecPushPull.shape[0]))
        
        for i in range(self._cmdHistory.shape[1]):
            j=2*i
            matrixToApply.T[j]= matrixWithAmp.T[i]* vecPushPull[0]
            matrixToApply.T[j+1]= matrixWithAmp.T[i]* vecPushPull[1]
        
        return matrixToApply
        
    
    def save(self):
        storeInFolder= CmdHistory._storageFolder()
        save= trackingNumberFolder.TtFolder(storeInFolder)
        dove, tt= save._createFolderToStoreMeasurements()
        
        fitsFileName= os.path.join(dove, 'info.fits')
        header= pyfits.Header()
        header['NPUSHPUL']= self._nPushPull
        header['WHO']= self._who
        pyfits.writeto(fitsFileName, self._modeVector, header)
        pyfits.append(fitsFileName, self._indexingList, header)
        pyfits.append(fitsFileName, self._cmdMatrix, header)
        pyfits.append(fitsFileName, self._cmdHToApply, header)
        pyfits.append(fitsFileName, self._ampVect, header)
        
        return tt
        
    @staticmethod
    def load(device, tt):
        theObject= CmdHistory(device)
        theObject._tt= tt
        storeInFolder= CmdHistory._storageFolder()
        folder= os.path.join(storeInFolder, tt)
        additionalInfoFitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.getheader(additionalInfoFitsFileName)
        hduList= pyfits.open(additionalInfoFitsFileName)
        theObject._modeVector= hduList[0].data
        theObject._indexingList= hduList[1].data
        theObject._cmdMatrix= hduList[2].data
        theObject._cmdHToApply= hduList[3].data
        theObject._ampVect= hduList[4].data
        theObject._who= header['WHO']
        try:
            theObject._nPushPull= header['NPUSHPUL']
        except KeyError:
            theObject._nPushPull= 1
            
        return theObject
        
        
        