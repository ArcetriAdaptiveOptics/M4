'''
@author: cs
'''
import os
import pyfits
import numpy as np
from m4.type.modalAmplitude import ModalAmplitude
from m4.type.modalBase import ModalBase
from m4.type.modesVector import ModesVector
from m4.ground.configuration import Configuration

def readObjectFitsFileName(amplitude_fits_file_name,
                           mode_vector_fits_file_name,
                           cmd_matrix_fits_file_name):
    '''
    return:
        Amplitude= vettore con l'ampiezza dei modi (numpy.array([]))
        modesVector= vettore dell'indice dei modi o degli attuatori
                    da applicare (numpy.array([]))
        cmd_matrix= matrice dei comandi dei modi
                    (nActs x nModes)
                     matrice diagonale nel caso di comandi zonali
    '''
    ma = ModalAmplitude.loadFromFits(amplitude_fits_file_name)
    amplitude = ma.getModalAmplitude()

    mv = ModesVector.loadFromFits(mode_vector_fits_file_name)
    mode_vector = mv.getModesVector()

    mb = ModalBase.loadFromFits(cmd_matrix_fits_file_name)
    cmd_matrix = mb.getModalBase()
    return amplitude, mode_vector, cmd_matrix

def readImageFromFitsFileName(fits_file_path):
    file_name = os.path.join(Configuration.TEST_IMAGE_ROOT_FOLDER,
                             fits_file_path)
    hduList = pyfits.open(file_name)
    ima = hduList[0].data
    immagine = np.ma.masked_array(ima[0], mask=np.invert(ima[1].astype(bool)))
    return immagine

def readDataFromFileFits(fits_file_path):
    file_name = os.path.join(Configuration.TEST_IMAGE_ROOT_FOLDER,
                             fits_file_path)
    hduList = pyfits.open(file_name)
    data = hduList[0].data
    return data
