'''
@author: cs
'''

import os

class Configuration(object):
        M4COORDINATE_ROOT_FOLDER= '/Users/rm/Desktop/Arcetri/M4/ActuatorCoordinates.fits'
        
        CALIBRATION_ROOT_FOLDER= "/Users/rm/Desktop/Arcetri/M4/ProvaCodice"
        
        FLATTENING_ROOT_FOLDER= os.path.join(Configuration.CALIBRATION_ROOT_FOLDER,
                                       "IFFunctions")