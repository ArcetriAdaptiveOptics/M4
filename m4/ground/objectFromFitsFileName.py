'''
@author: cs
'''

from m4.type.modalAmplitude import ModalAmplitude
from m4.type.modalBase import ModalBase
from m4.type.modesVector import ModesVector
from m4.ground.configuration import Configuration
import os
import pyfits
import numpy as np

def readObjectFitsFileName(amplitudeFitsFileName, modeVectorFitsFileName, cmdMatrixFitsFileName):
        '''
        return:
            Amplitude= vettore con l'ampiezza dei modi (numpy.array([]))
            modesVector= vettore dell'indice dei modi o degli attuatori
                        da applicare (numpy.array([]))
            cmdMatrix= matrice dei comandi dei modi 
                         (nActs x nModes)
                         matrice diagonale nel caso di comandi zonali
        '''
        ma= ModalAmplitude.loadFromFits(amplitudeFitsFileName)
        amplitude= ma.getModalAmplitude()
        
        mv= ModesVector.loadFromFits(modeVectorFitsFileName)
        modeVector= mv.getModesVector()
        
        mb= ModalBase.loadFromFits(cmdMatrixFitsFileName)
        cmdMatrix= mb.getModalBase()
        return amplitude, modeVector, cmdMatrix
    
    
def readImageFromFitsFileName(fitsFilePath):
    fileName= os.path.join(Configuration.TEST_IMAGE_ROOT_FOLDER, fitsFilePath)
    hduList= pyfits.open(fileName)
    ima= hduList[0].data
    immagine= np.ma.masked_array(ima[0], mask= np.invert(ima[1].astype(bool)))
    return immagine

def readDataFromFileFits(fitsFilePath):
    fileName= os.path.join(Configuration.TEST_IMAGE_ROOT_FOLDER, fitsFilePath)
    hduList= pyfits.open(fileName)
    data= hduList[0].data
    return data