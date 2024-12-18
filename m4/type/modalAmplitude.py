'''
Autors
  - C. Selmi: written in 2019
'''

import os
import numpy as np
from astropy.io import fits as pyfits
import h5py
from m4.configuration import config_folder_names as fold_name


class ModalAmplitude():
    '''
    Class to create and manage the modal amplitude

    HOW TO USE IT::

        from m4.type.modalAmplitude import ModalAmplitude
        ma = ModalAmplitude()
    '''

    def __init__(self):
        """The constructor """
        self._modalAmplitude = None
        self._fitsfilename = None
        self.tag = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save data"""
        return fold_name.MODALAMPLITUDE_ROOT_FOLDER

    def getModalAmplitude(self):
        '''
        Returns
        -------
        modalAmplitude: numpy array
                        vector of modal amplitude
        '''
        return self._modalAmplitude

    def getTag(self):
        '''
        Returns
        -------
        tag: string
            modal ampplitude tag
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


    def saveAsFits(self, tag, modal_amplitude):
        ''' Save the data in fits format

        Parameters
        ----------
            tag: string
                file name to save
            modal_amplitude: numpy array
                            vector of amplitude
        '''
        self._tag = tag
        store_in_folder = ModalAmplitude._storageFolder()
        filename = tag + '.fits'
        fits_file_name = os.path.join(store_in_folder, filename)
        pyfits.writeto(fits_file_name, modal_amplitude)

    def saveAsH5(self, tag, modal_amplitude):
        ''' Save the data in h5 format

        Parameters
        ----------
            tag: string
                file name to save
            modal_amplitude: numpy array
                            vector of amplitude
        '''
        store_in_folder = ModalAmplitude._storageFolder()
        filename = tag + '.h5'
        hf = h5py.File(os.path.join(store_in_folder, filename), 'w')
        hf.create_dataset('dataset_1', data=modal_amplitude)
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
                    ModalAmplitude class object
        """
        theObject = ModalAmplitude()
        store_in_folder = ModalAmplitude._storageFolder()
        all_fits_file_name = os.path.join(store_in_folder, fits_file_name)
        hduList = pyfits.open(all_fits_file_name)
        theObject._modalAmplitude = hduList[0].data
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
                    ModalAmplitude class object
        """
        store_in_folder = ModalAmplitude._storageFolder()
        hf = h5py.File(os.path.join(store_in_folder, filename), 'r')
        hf.keys()
        data = hf.get('dataset_1')
        modal_amplitude = np.array(data)
        hf.close()
        theObject = ModalAmplitude()
        theObject._modalAmplitude = modal_amplitude
        theObject._fitsfilename = filename
        return theObject
