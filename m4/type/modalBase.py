'''
@author: cs
'''

import os
from astropy.io import fits as pyfits
import h5py
import numpy as np
from m4.configuration import config_folder_names as fold_name


class ModalBase():
    '''
    Class to create and manage the modal base

    HOW TO USE IT::

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
        return fold_name.MODALBASE_ROOT_FOLDER

    def getModalBase(self):
        '''
        Returns
        -------
        modalBase: numpy array
                    vector of modal base
        '''
        return self._modalBase

    def getTag(self):
        '''
        Returns
        -------
        tag: string
            modal base tag
        '''
        return self._tag

    def getFitsFileName(self):
        '''
        Returns
        -------
        fitsfilename: string
                    path fits file name
        '''
        return self._fitsfilename

    def saveAsFits(self, tag, modal_base):
        ''' Save the data in fits format

        Parameters
        ----------
            tag: string
                file name to save
            modal_base: numpy array [nActs x nModes]
                            command matrix
        '''
        self._tag = tag
        store_in_folder = ModalBase._storageFolder()
        filename = tag + '.fits'
        fits_file_name = os.path.join(store_in_folder, filename)
        pyfits.writeto(fits_file_name, modal_base)

    def saveAsH5(self, tag, modal_base):
        ''' Save the data in h5 format

        Parameters
        ----------
            tag: string
                file name to save
            modal_bace: numpy array
                        vector of modal base
        '''
        store_in_folder = ModalBase._storageFolder()
        filename = tag + '.h5'
        hf = h5py.File(os.path.join(store_in_folder, filename), 'w')
        hf.create_dataset('dataset_1', data=modal_base)
        hf.close()

    @staticmethod
    def loadFromFits(fits_file_name):
        """ Creates the object

        Parameters
        ----------
        fits_file_name : string
                        modal amplitude path

        Returns
        -------
        theObject: object
                    ModalBase class object
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

        Parameters
        ----------
        filename : string
                 modal amplitude path

        Returns
        -------
        theObject: object
                    ModalBase class object
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
    