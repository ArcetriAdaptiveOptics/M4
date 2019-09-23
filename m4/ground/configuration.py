'''
@author: cs
'''

import os
import numpy as np

class Configuration(object):
    #ROOT
        M4COORDINATE_ROOT_FOLDER= '/Users/rm/Desktop/Arcetri/M4/ActuatorCoordinates.fits'
        
        CALIBRATION_ROOT_FOLDER= "/Users/rm/Desktop/Arcetri/M4/ProvaCodice"
        
        
    #Interf  
        Lambda= 2
    #Ott
        #ParabolaPupilXYRadius= np.array([128, 128, 128]) #centerX, centerY, radius
        ParabolaPupilXYRadius= np.array([256, 256, 256])
    #DM
        n_Seg= 6
        nActsTot= 5352
        nActSeg= 892