'''
@author: cs
'''

import os
from astropy.io import fits as pyfits
import h5py
import numpy as np
from m4.configuration import config_folder_names as fold_name


class ModesVector(object):
    '''
    Class to create and manage the mode vector

    HOW TO USE IT::

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
        return fold_name.MODESVECTOR_ROOT_FOLDER

    def getModesVector(self):
        '''
        Returns
        -------
        modesVector: numpy array
                    vector of modes
        '''
        return self._modesVector

    def getTag(self):
        '''
        Returns
        -------
        tag: string
            mode vector tag
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


    def saveAsFits(self, tag, modes_vector):
        ''' Save the data in fits format

        Parameters
        ----------
            tag: string
                file name to save
            modes_vector: numpy array
                            vector of selected modes
        '''
        self._tag = tag
        store_in_folder = ModesVector._storageFolder()
        filename = tag + '.fits'
        fits_file_name = os.path.join(store_in_folder, filename)
        pyfits.writeto(fits_file_name, modes_vector)

    def saveAsH5(self, tag, modes_vector):
        ''' Save the data in h5 format

        Parameters
        ----------
            tag: string
                file name to save
            modes_vector: numpy array
                        vector of selectred modes
        '''
        store_in_folder = ModesVector._storageFolder()
        filename = tag + '.h5'
        hf = h5py.File(os.path.join(store_in_folder, filename), 'w')
        hf.create_dataset('dataset_1', data=modes_vector)
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
                    ModaesVector class object
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

        Parameters
        ----------
        filename : string
                 modal amplitude path

        Returns
        -------
        theObject: object
                    ModesVector class object
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
    