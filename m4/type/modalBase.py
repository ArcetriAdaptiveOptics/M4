'''
@author: cs
'''

import os
from astropy.io import fits as pyfits
import h5py
import numpy as np
from m4.ground.configuration import Configuration


class ModalBase():
    '''
    Class to create and manage the modal base

    HOW TO USE IT:
    from m4.type.modalBase import ModalBase
    mb = ModalBase()
    '''

    def __init__(self):
        """The constructor """
        self._modalBase = None
        self._fitsfilename = None
        self._tag = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save data"""
        return os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                            "ModalBase")

    def getModalBase(self):
        return self._modalBase

    def getTag(self):
        return self._tag

    def getFitsFileName(self):
        return self._fitsfilename

    def saveAsFits(self, tag, modal_base):
        self._tag = tag
        ''' Save the data in h5 format
        args:
            tag (string) = file name to save
            modal_base = command matrix
                                (nActs x nModes)
        '''
        store_in_folder = ModalBase._storageFolder()
        filename = tag + '.fits'
        fits_file_name = os.path.join(store_in_folder, filename)
        pyfits.writeto(fits_file_name, modal_base)

    def saveAsH5(self, tag, modal_base):
        ''' Save the data in h5 format
        args:
            tag (string) = file name to save
            modes_base = command matrix
        '''
        store_in_folder = ModalBase._storageFolder()
        filename = tag + '.h5'
        hf = h5py.File(os.path.join(store_in_folder, filename), 'w')
        hf.create_dataset('dataset_1', data=modal_base)
        hf.close()

    @staticmethod
    def loadFromFits(fits_file_name):
        """ Creates the object

            Args:
                fits_file_name = modal base path

            Returns:
                    theObject = ModalBase class object
        """
        theObject = ModalBase()
        store_in_folder = ModalBase._storageFolder()
        all_fits_file_name = os.path.join(store_in_folder, fits_file_name)
        hduList = pyfits.open(all_fits_file_name)
        theObject._modalBase = hduList[0].data
        theObject._fitsfilename = fits_file_name
        return theObject

    @staticmethod
    def loadFromH5(filename):
        """ Creates the object

            Args:
                filename = modal base path

            Returns:
                    theObject = ModalBase class object
        """
        theObject = ModalBase()
        theObject._fitsfilename = filename
        store_in_folder = ModalBase._storageFolder()
        hf = h5py.File(os.path.join(store_in_folder, filename), 'r')
        hf.keys()
        data = hf.get('dataset_1')
        theObject._modalBase = np.array(data)
        hf.close()
        return theObject
    