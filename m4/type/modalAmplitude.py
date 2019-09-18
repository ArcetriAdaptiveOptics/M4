'''
@author: cs
'''


from m4.ground.configuration import Configuration
import os   
import pyfits
import h5py
import numpy as np
    
    
    
class ModalAmplitude():
    
    def __init__(self):
        self._modalAmplitude= None
        self._fitsfilename= None
    
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "ModalAmplitude")
    
    def getModalAmplitude(self):
        return self._modalAmplitude
    
    
    def saveAsFits(self, tag, modalAmplitude):
        '''
            tag (stringa)= nome del file da salvare
            modalAmplitude (array)= vettore delle ampiezze
        '''
        storeInFolder= ModalAmplitude._storageFolder()
        filename= tag + '.fits'
        fitsFileName= os.path.join(storeInFolder, filename)
        pyfits.writeto(fitsFileName, modalAmplitude)
        
    def saveAsH5(self, tag, modalAmplitude):
        storeInFolder= ModalAmplitude._storageFolder()
        filename= tag + '.h5'
        hf = h5py.File(os.path.join(storeInFolder,filename), 'w')
        hf.create_dataset('dataset_1', data=modalAmplitude)
        hf.close()
    
    @staticmethod 
    def loadFromFits(fitsfilename):
        theObject= ModalAmplitude()
        storeInFolder= ModalAmplitude._storageFolder()
        allFitsFileName= os.path.join(storeInFolder, fitsfilename)
        hduList= pyfits.open(allFitsFileName)
        theObject._modalAmplitude= hduList[0].data
        theObject._fitsfilename= fitsfilename
        return theObject
    
    @staticmethod 
    def loadFromH5(filename):
        storeInFolder= ModalAmplitude._storageFolder()
        hf = h5py.File(os.path.join(storeInFolder,filename), 'r')
        hf.keys()
        data= hf.get('dataset_1')
        modalAmplitude= np.array(data)
        hf.close()
        theObject= ModalAmplitude()
        theObject._modalAmplitude= modalAmplitude
        theObject._fitsfilename= filename
        return theObject
        
        