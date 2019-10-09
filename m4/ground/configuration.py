'''
@author: cs
'''

import os
import numpy as np

class Configuration(object):
    #ROOT
        LOG_ROOT_FOLDER= '/Users/rm/Desktop/Arcetri/M4/ProvaCodice'
        
        M4COORDINATE_ROOT_FOLDER= '/Users/rm/Desktop/Arcetri/M4/ActuatorCoordinates.fits'
        
        CALIBRATION_ROOT_FOLDER= "/Users/rm/Desktop/Arcetri/M4/ProvaCodice"
        
        TEST_IMAGE_ROOT_FOLDER= "/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova"
        
        
    #Interf  
        Lambda= 2
    #Ott
        #ParabolaPupilXYRadius= np.array([128, 128, 128]) #centerX, centerY, radius
        ParabolaPupilXYRadius= np.array([256, 256, 256])
        ParabolDoF= 3
        RMDof= 2
    #DM
        n_Seg= 6
        nActsTot= 5352
        nActSeg= 892
        M4Dof= 6