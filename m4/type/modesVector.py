'''
@author: cs
'''

from m4.ground.configuration import Configuration
import os   
import pyfits
import h5py
import numpy as np


class ModesVector(object):

    def __init__(self):
        self._modesVector= None 
        self._fitsfilename= None
        self.tag= None 
        
    @staticmethod
    def _storageFolder():
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "ModesVector")
        
    def getModesVector(self):
        return self._modesVector
    
    def getTag(self):
        return self._tag
    
    def getFitsFileName(self):
        return self._fitsfilename
    
    
    def saveAsFits(self, tag, modesVector):
        self._tag= tag
        '''
            tag (stringa)= nome del file da salvare
            modesVector= vettore dei modi scelti
        '''
        storeInFolder= ModesVector._storageFolder()
        filename= tag + '.fits'
        fitsFileName= os.path.join(storeInFolder, filename)
        pyfits.writeto(fitsFileName, modesVector)
        
    def saveAsH5(self, tag, modesVector):
        storeInFolder= ModesVector._storageFolder()
        filename= tag + '.h5'
        hf = h5py.File(os.path.join(storeInFolder,filename), 'w')
        hf.create_dataset('dataset_1', data=modesVector)
        hf.close()
    
    @staticmethod 
    def loadFromFits(fitsfilename):
        theObject= ModesVector()
        storeInFolder= ModesVector._storageFolder()
        allFitsFileName= os.path.join(storeInFolder, fitsfilename)
        hduList= pyfits.open(allFitsFileName)
        theObject._modesVector= hduList[0].data
        theObject._fitsfilename= fitsfilename
        return theObject
    
    @staticmethod 
    def loadFromH5(filename):
        theObject= ModesVector()
        theObject._fitsfilename= filename
        storeInFolder= ModesVector._storageFolder()
        hf = h5py.File(os.path.join(storeInFolder,filename), 'r')
        hf.keys()
        data= hf.get('dataset_1')
        theObject._modesVector= np.array(data)
        hf.close()
        return theObject