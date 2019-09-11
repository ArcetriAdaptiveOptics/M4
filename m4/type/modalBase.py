'''
@author: cs
'''

from m4.ground.configuration import Configuration
import os 
import pyfits


class ModalBase():
    
    def __init__(self):
        self._modalBase= None 
        self._fitsfilename= None
    
    def getModalBase(self):
        return self.modalBase
    
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "ModalBase")
        
    def save(self, tag, modalBase):
        storeInFolder= ModalBase._storageFolder()
        filename= tag + '.fits'
        fitsFileName= os.path.join(storeInFolder, filename)
        pyfits.writeto(fitsFileName, modalBase)
    
    @staticmethod 
    def load(fitsfilename):
        theObject= ModalBase()
        storeInFolder= ModalBase._storageFolder()
        allFitsFileName= os.path.join(storeInFolder, fitsfilename)
        hduList= pyfits.open(allFitsFileName)
        theObject._modalBase= hduList[0].data
        theObject._fitsfilename= fitsfilename
        
        return theObject