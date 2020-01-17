'''
@author: cs
'''
import os
from astropy.io import fits as pyfits
import numpy as np
from m4.type.modalAmplitude import ModalAmplitude
from m4.type.modalBase import ModalBase
from m4.type.modesVector import ModesVector
from m4.ground.configuration import Configuration

def readObjectFitsFileName(amplitude_fits_file_name,
                           mode_vector_fits_file_name,
                           cmd_matrix_fits_file_name):
    '''
    args:
        amplitude_fits_file_name = vector with mode amplitude fits file name
        mode_vector_fits_file_name = mode or actuator index vector
                                    to be applied fits file name
        cmd_matrix_fits_file_name = matrix of mode commands fits file name

    returns:
        amplitude = vector with mode amplitude (numpy.array([]))
        modesVector = mode or actuator index vector
                    to be applied (numpy.array([]))
        cmd_matrix = matrix of mode commands
                    (nActs x nModes)
                     diagonal matrix in case of zonal commands
    '''
    ma = ModalAmplitude.loadFromFits(amplitude_fits_file_name)
    amplitude = ma.getModalAmplitude()

    mv = ModesVector.loadFromFits(mode_vector_fits_file_name)
    mode_vector = mv.getModesVector()

    mb = ModalBase.loadFromFits(cmd_matrix_fits_file_name)
    cmd_matrix = mb.getModalBase()
    return amplitude, mode_vector, cmd_matrix

def readImageFromFitsFileName(fits_file_path):
    """
    args:
        fits_file_path = file path of the image.fits to read

    returns:
            immagine = masked array of the image
    """
    file_name = os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                             fits_file_path)
    hduList = pyfits.open(file_name)
    ima = hduList[0].data
    immagine = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
    return immagine

def readDataFromFileFits(fits_file_path):
    """
    args:
        fits_file_path = file path of the data to read

    returns:
            data = data included in file path
    """
    file_name = os.path.join(Configuration.OPD_DATA_FOLDER,
                             fits_file_path)
    hduList = pyfits.open(file_name)
    data = hduList[0].data
    return data

def readImageFromRunaIFFs(fits_file_path):
    file_name = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER,
                             fits_file_path)
    hduList = pyfits.open(file_name)
    cube = hduList[0].data
    immagine = np.ma.masked_array(cube[0,:,:], mask=np.invert(cube[1,:,:].astype(bool)))
    return immagine
