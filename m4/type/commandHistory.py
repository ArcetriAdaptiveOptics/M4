'''
@author: cs
'''

from m4.utils import saveIFFInfo
from m4.utils.configuration import Configuration
import numpy as np
import os 
import pyfits


class CmdHistory(object):
    
    def __init__(self, device):
        self._device= device[0]
        self._who= device[1]
        self._nActs = self._device.nActs()
        self._indexing= None 
        self._nPushPull= None
        self._cmdMatrix= None 
        
    
    def createShuffleCmdHistory(self, indexing, nPushPull, cmdMatrix):
        self._indexing= indexing
        self._nPushPull= nPushPull
        self._cmdMatrix= cmdMatrix
        
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
                
        self._matrixToApply= matrixToApply
        self._indexingList= np.array(indexingList)
                
        return matrixToApply, indexingList
    
    def createCmdHistory(self):
        pass
    
    
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                   "CommandHistory")


    def saveCmdHistory(self):
        storeInFolder= CmdHistory._storageFolder()
        save= saveIFFInfo.SaveAdditionalInfo(storeInFolder)
        dove= save._createFolderToStoreMeasurements()
        fitsFileName= os.path.join(dove, 'info.fits')
        header= pyfits.Header()
        header['NPUSHPUL']= self._nPushPull
        header['WHO']= self._who
        pyfits.writeto(fitsFileName, self._indexing, header)
        pyfits.append(fitsFileName, self._indexingList, header)
        pyfits.append(fitsFileName, self._cmdMatrix, header)

    @staticmethod
    def loadCmdHistory(device, tt):
        theObject= CmdHistory(device)
        theObject._tt= tt
        storeInFolder= CmdHistory._storageFolder()
        folder= os.path.join(storeInFolder, tt)
        additionalInfoFitsFileName= os.path.join(folder, 'info.fits')
        header= pyfits.getheader(additionalInfoFitsFileName)
        hduList= pyfits.open(additionalInfoFitsFileName)
        theObject._indexing= hduList[0].data
        theObject._indexingList= hduList[1].data
        theObject._cmdMatrix= hduList[2].data
        theObject._who= header['WHO']
        try:
            theObject._nPushPull= header['NPUSHPUL']
        except KeyError:
            theObject._nPushPull= 1
            
        return theObject
        
        
        