"""
Autors
  - C. Selmi:  written in 2019

Function for reading file fits::

    from m4.ground import object_from_fits_file_name as obj
    numpy_array_data = obj.readFits_object(fits_file_path)
    or
    numpy_masked_array_data = obj.readFits_maskedImage(fits_file_path)
    or
    amplitude, mode_vector, cmd_matrix = obj.readObjectFitsFileName(
                                'ampName.fits', 'mvec.fits', 'cmdMatrix.fits')
"""
import os
from astropy.io import fits as pyfits
import numpy as np
from m4.type.modalAmplitude import ModalAmplitude
from m4.type.modalBase import ModalBase
from m4.type.modesVector import ModesVector
from m4.configuration.config import *

### Generiche
def readFits_object(fits_file_path):
    '''
    Parameters
    ----------
    fits_file_path: string
                    fits file path of numpy array to read

    Returns
    -------
            object: numpy array
    '''
    hduList = pyfits.open(fits_file_path)
    obj = hduList[0].data
    return obj

def readFits_maskedImage(fits_file_path):
    '''
    Parameters
    ----------
    fits_file_path: string
                    fits file path of masked array to read

    Returns
    -------
    object: numpy masked array
    '''
    hduList = pyfits.open(fits_file_path)
    ima = hduList[0].data
    immagine = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
    return immagine
###

def _readImageFromFitsFileName(fits_file_path):
    """
    Parameters
    ----------
        fits_file_path : string
                        file path of the image.fits to read

    Returns
    -------
            immagine: numpy masked array
                        masked array of the image
    """
    file_name = os.path.join(path_name.CALIBRATION_ROOT_FOLDER,
                             fits_file_path)
    hduList = pyfits.open(file_name)
    ima = hduList[0].data
    immagine = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
    return immagine

def readObjectFitsFileName(amplitude_fits_file_name,
                           mode_vector_fits_file_name,
                           cmd_matrix_fits_file_name):
    '''
    Parameters
    ----------
        amplitude_fits_file_name: string
                                 vector with mode amplitude fits file name
        mode_vector_fits_file_name: string
                                    mode or actuator index vector
                                    to be applied fits file name
        cmd_matrix_fits_file_name: string
                                matrix of mode commands fits file name

    Returns
    -------
        amplitude: numpy array
                vector with mode amplitude
        modesVector: numpy array
                    mode or actuator index vector to be applied
        cmd_matrix: numpy array [nActs x nModes]
                    matrix of mode commands
                    diagonal matrix in case of zonal commands
    '''
    ma = ModalAmplitude.loadFromFits(amplitude_fits_file_name)
    amplitude = ma.getModalAmplitude()

    mv = ModesVector.loadFromFits(mode_vector_fits_file_name)
    mode_vector = mv.getModesVector()

    mb = ModalBase.loadFromFits(cmd_matrix_fits_file_name)
    cmd_matrix = mb.getModalBase()
    return amplitude, mode_vector, cmd_matrix

def _readDataFromFileFits(fits_file_path):
    """
    Parameters
    ----------
        fits_file_path: string
                    file path of the data to read

    Returns
    -------
            data: numpy array
                    data included in file path
    """
    file_name = os.path.join(path_name.OPT_DATA_FOLDER,
                             fits_file_path)
    hduList = pyfits.open(file_name)
    data = hduList[0].data
    return data

def _readImageFromRunaIFFs(fits_file_path):
    file_name = os.path.join(fold_name.IFFUNCTIONS_ROOT_FOLDER,
                             fits_file_path)
    hduList = pyfits.open(file_name)
    cube = hduList[0].data
    immagine = np.ma.masked_array(cube[0, :, :], mask=np.invert(cube[1, :, :].astype(bool)))
    return immagine
