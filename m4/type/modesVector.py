'''
@author: cs
'''

import os
import pyfits
import h5py
import numpy as np
from m4.ground.configuration import Configuration


class ModesVector(object):
    '''
    Classe per creare e gestire il vettore dei modi
    '''

    def __init__(self):
        self._modesVector = None
        self._fitsfilename = None
        self.tag = None

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


    def saveAsFits(self, tag, modes_vector):
        self._tag = tag
        '''
            tag (stringa)= nome del file da salvare
            modes_vector= vettore dei modi scelti
        '''
        store_in_folder = ModesVector._storageFolder()
        filename = tag + '.fits'
        fits_file_name = os.path.join(store_in_folder, filename)
        pyfits.writeto(fits_file_name, modes_vector)

    def saveAsH5(self, tag, modes_vector):
        store_in_folder = ModesVector._storageFolder()
        filename = tag + '.h5'
        hf = h5py.File(os.path.join(store_in_folder, filename), 'w')
        hf.create_dataset('dataset_1', data=modes_vector)
        hf.close()

    @staticmethod
    def loadFromFits(fits_file_name):
        theObject = ModesVector()
        store_in_folder = ModesVector._storageFolder()
        all_fits_file_name = os.path.join(store_in_folder, fits_file_name)
        hduList = pyfits.open(all_fits_file_name)
        theObject._modesVector = hduList[0].data
        theObject._fitsfilename = fits_file_name
        return theObject

    @staticmethod
    def loadFromH5(filename):
        theObject = ModesVector()
        theObject._fitsfilename = filename
        store_in_folder = ModesVector._storageFolder()
        hf = h5py.File(os.path.join(store_in_folder, filename), 'r')
        hf.keys()
        data = hf.get('dataset_1')
        theObject._modesVector = np.array(data)
        hf.close()
        return theObject
    