'''
@author: cs
'''

from m4.ground import trackingNumberFolder
from m4.ground.configuration import Configuration
import numpy as np
from m4.ground import logger
import os 
import copy
import pyfits
import h5py


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
    
    def getIndexingList(self):
        return self._indexingList
        
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
         return:
             matrixToApply: command history (nAct, nModes x nPushPusll x 2)
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
        tt= self.saveAsFits()
        logger.log('Creazione della commandHistoryMatrix', 'ordinata', tt)
        print(tt)
        
        return matrixToApply, tt
    
    
    def shuffleCommandHistoryMaker(self, modeVector, ampVector, 
                                        cmdMatrix, nPushPull):
        self._ampVect= ampVector
        cmdHistory, indexingList= self._shuffleCmdHistory(
                                        modeVector, nPushPull, cmdMatrix)
        zipped= self._zippedAmplitude(ampVector)
        matrixToApply= self._cmdHistoryToApply(zipped)
        self._cmdHToApply=matrixToApply
        tt= self.saveAsFits()
        logger.log('Creazione della commandHistoryMatrix', 'casuale', tt)
        print(tt)
        
        return matrixToApply, tt
        
    
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
        indList=[]
        for i in range(nPushPull):
            indList.append(modeVector)
        self._indexingList= np.array(indList)  
        #self._indexingList= np.tile(modeVector, nPushPull)
        
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
        
    
    def saveAsFits(self):
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
    
    def saveAsH5(self):
        storeInFolder= CmdHistory._storageFolder()
        save= trackingNumberFolder.TtFolder(storeInFolder)
        dove, tt= save._createFolderToStoreMeasurements()
        
        fitsFileName= os.path.join(dove, 'info.h5')
        hf = h5py.File(fitsFileName, 'w')
        hf.create_dataset('dataset_1', data=self._modeVector)
        hf.create_dataset('dataset_2', data=self._indexingList)
        hf.create_dataset('dataset_3', data=self._cmdHToApply)
        hf.create_dataset('dataset_4', data=self._ampVect)
        hf.attrs['NPUSHPUL']= self._nPushPull
        hf.attrs['WHO']= self._who
        hf.close()
        return tt
        
    @staticmethod
    def loadFromFits(device, tt):
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
    
    @staticmethod
    def loadFromH5(device, tt):
        theObject= CmdHistory(device)
        theObject._tt= tt
        storeInFolder= CmdHistory._storageFolder()
        folder= os.path.join(storeInFolder, tt)
        fileName= os.path.join(folder, 'info.h5')
        hf = h5py.File(fileName, 'r')
        hf.keys()
        data1= hf.get('dataset_1')
        data2= hf.get('dataset_2')
        data3= hf.get('dataset_3')
        data4= hf.get('dataset_4')
        theObject._nPushPull= hf.attrs['NPUSHPUL']
        theObject._who= hf.attrs['WHO']
        theObject._modeVector= np.array(data1)
        theObject._indexingList= np.array(data2)
        theObject._cmdHToApply= np.array(data3)
        theObject._ampVect= np.array(data4)
        hf.close()
        return theObject
        
        
        
        
        