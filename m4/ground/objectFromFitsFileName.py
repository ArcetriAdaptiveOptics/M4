'''
@author: cs
'''

from m4.type.modalAmplitude import ModalAmplitude
from m4.type.modalBase import ModalBase
from m4.type.modesVector import ModesVector
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
    
    
def readImageFromFitsFileName(fitsFileName):
    fitsRoot= "/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova"
    fileName= os.path.join(fitsRoot, fitsFileName)
    hduList= pyfits.open(fileName)
    ima= hduList[0].data
    immagine= np.ma.masked_array(ima[0], mask= ima[1])
    return immagine