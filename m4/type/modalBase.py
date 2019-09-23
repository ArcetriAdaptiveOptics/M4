'''
@author: cs
'''

from m4.ground.configuration import Configuration
import os 
import pyfits
import h5py
import numpy as np 


class ModalBase():
    
    def __init__(self):
        self._modalBase= None 
        self._fitsfilename= None
        self._tag= None
    
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "ModalBase")
        
    def getModalBase(self):
        return self._modalBase

    def getTag(self):
        return self._tag
    
    def getFitsFileName(self):
        return self._fitsfilename
        
        
    def saveAsFits(self, tag, modalBase):
        self._tag= tag
        '''
            tag (stringa)= nome del file da salvare
            modalBase= matrice dei comandi 
                                (nActs x nModes)
        '''
        storeInFolder= ModalBase._storageFolder()
        filename= tag + '.fits'
        fitsFileName= os.path.join(storeInFolder, filename)
        pyfits.writeto(fitsFileName, modalBase)
        
    def saveAsH5(self, tag, modalBase):
        storeInFolder= ModalBase._storageFolder()
        filename= tag + '.h5'
        hf = h5py.File(os.path.join(storeInFolder,filename), 'w')
        hf.create_dataset('dataset_1', data=modalBase)
        hf.close()
    
    @staticmethod 
    def loadFromFits(fitsfilename):
        theObject= ModalBase()
        storeInFolder= ModalBase._storageFolder()
        allFitsFileName= os.path.join(storeInFolder, fitsfilename)
        hduList= pyfits.open(allFitsFileName)
        theObject._modalBase= hduList[0].data
        theObject._fitsfilename= fitsfilename
        return theObject
    
    @staticmethod 
    def loadFromH5(filename):
        theObject= ModalBase()
        theObject._fitsfilename= filename
        storeInFolder= ModalBase._storageFolder()
        hf = h5py.File(os.path.join(storeInFolder,filename), 'r')
        hf.keys()
        data= hf.get('dataset_1')
        theObject._modalBase= np.array(data)
        hf.close()
        return theObject