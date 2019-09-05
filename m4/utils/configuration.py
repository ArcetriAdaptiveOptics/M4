'''
@author: cs
'''

import os
import numpy as np

class Configuration(object):
    #ROOT
        M4COORDINATE_ROOT_FOLDER= '/Users/rm/Desktop/Arcetri/M4/ActuatorCoordinates.fits'
        
        CALIBRATION_ROOT_FOLDER= "/Users/rm/Desktop/Arcetri/M4/ProvaCodice"
        
        FLATTENING_ROOT_FOLDER= os.path.join(CALIBRATION_ROOT_FOLDER, "IFFunctions")
        
    #Interf  
        Lambda= 633
    #Ott
        ParabolaPupilXYRadius= np.array([0, 0, 1]) #centerX, centerY, radius
    #DM
        n_Seg= 6
        nActsTot= 5352
        nActSeg= 892