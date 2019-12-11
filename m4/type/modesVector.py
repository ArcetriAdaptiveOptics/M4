'''
@author: cs
'''

import os
from astropy.io import fits as pyfits
import h5py
import numpy as np
from m4.ground.configuration import Configuration


class ModesVector(object):
    '''
    Class to create and manage the mode vector

    HOW TO USE IT:
    from m4.type.modesVector import ModesVestor
    mv = ModesVector()
    '''

    def __init__(self):
        """The constructor """
        self._modesVector = None
        self._fitsfilename = None
        self.tag = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save data"""
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "ModesVector")

    def getModesVector(self):
        return self._modesVector

    def getTag(self):
        return self._tag

    def getFitsFileName(self):
        return self._fitsfilename


    def saveAsFits(self, tag, modes_vector):
        ''' Save the data in fits format
        args:
            tag (string) = file name to save
            modes_vector = vector of selected modes
        '''
        self._tag = tag
        store_in_folder = ModesVector._storageFolder()
        filename = tag + '.fits'
        fits_file_name = os.path.join(store_in_folder, filename)
        pyfits.writeto(fits_file_name, modes_vector)

    def saveAsH5(self, tag, modes_vector):
        ''' Save the data in h5 format
        args:
            tag (string) = file name to save
            modes_vector = vector of selected modes
        '''
        store_in_folder = ModesVector._storageFolder()
        filename = tag + '.h5'
        hf = h5py.File(os.path.join(store_in_folder, filename), 'w')
        hf.create_dataset('dataset_1', data=modes_vector)
        hf.close()

    @staticmethod
    def loadFromFits(fits_file_name):
        """ Creates the object

            Args:
                filename = mode vector path

            Returns:
                    theObject = ModalVector class object
        """
        theObject = ModesVector()
        store_in_folder = ModesVector._storageFolder()
        all_fits_file_name = os.path.join(store_in_folder, fits_file_name)
        hduList = pyfits.open(all_fits_file_name)
        theObject._modesVector = hduList[0].data
        theObject._fitsfilename = fits_file_name
        return theObject

    @staticmethod
    def loadFromH5(filename):
        """ Creates the object

            Args:
                filename = mode vector path

            Returns:
                    theObject = ModalVector class object
        """
        theObject = ModesVector()
        theObject._fitsfilename = filename
        store_in_folder = ModesVector._storageFolder()
        hf = h5py.File(os.path.join(store_in_folder, filename), 'r')
        hf.keys()
        data = hf.get('dataset_1')
        theObject._modesVector = np.array(data)
        hf.close()
        return theObject
    