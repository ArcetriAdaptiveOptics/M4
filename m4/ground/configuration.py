'''
@author: cs
'''
import numpy as np


class Configuration():
    '''
    Costanti di configurazione del sistema
    '''
    #ROOT
    LOG_ROOT_FOLDER = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice'

    M4COORDINATE_ROOT_FOLDER = \
        '/Users/rm/Desktop/Arcetri/M4/ActuatorCoordinates.fits'

    CALIBRATION_ROOT_FOLDER = "/Users/rm/Desktop/Arcetri/M4/ProvaCodice"

    TEST_IMAGE_ROOT_FOLDER = \
        "/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova"


    #Interf
    LAMBA = 2
    #Ott
    #ParabolaPupilXYRadius= np.array([128, 128, 128]) #centerX, centerY, radius
    PARABOLA_PUPIL_XYRADIUS = np.array([256, 256, 256])
    PARABOLA_DOF = 3
    RM_DOF = 2
    #DM
    N_SEG = 6
    N_ACTS_TOT = 5352
    N_ACT_SEG = 892
    M4_DOF = 6
