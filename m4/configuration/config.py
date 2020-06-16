'''
@author: cs
'''
import os

class path_name():
    '''
    '''
    BASE_PATH = '/mnt/shared/'
    CONFIGURATION_ROOT_FOLDER = os.path.join(BASE_PATH, 'SYSCONFData')
    CALIBRATION_ROOT_FOLDER = os.path.join(BASE_PATH, 'M4Data')
    OPD_DATA_FOLDER = os.path.join(CALIBRATION_ROOT_FOLDER, 'OPTData')
    OUT_FOLDER = os.path.join(CALIBRATION_ROOT_FOLDER, 'Results')

class fold_name(path_name):
    '''
    '''
    LOG_ROOT_FOLDER = os.path.join(path_name.CALIBRATION_ROOT_FOLDER, 'LOGData/mylog')
    IFFUNCTIONS_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER, 'IFFunctions')
