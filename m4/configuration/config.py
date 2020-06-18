'''
@author: cs
'''
import os

mirror_conf = '20170430'
optical_conf= '20150730'
simulated = 1 

class path_name():
    '''
    '''
    #BASE_PATH = '/mnt/shared/'
    BASE_PATH = '/Users/rm/Desktop/Arcetri/M4/Data'
    CONFIGURATION_ROOT_FOLDER = os.path.join(BASE_PATH, 'SYSCONFData')
    CALIBRATION_ROOT_FOLDER = os.path.join(BASE_PATH, 'M4Data')
    OPD_DATA_FOLDER = os.path.join(CALIBRATION_ROOT_FOLDER, 'OPTData')
    OUT_FOLDER = os.path.join(CALIBRATION_ROOT_FOLDER, 'Results')
    MIRROR_FOLDER = os.path.join(BASE_PATH, 'MIRROR_System')
    OPTICAL_FOLDER = os.path.join(BASE_PATH, 'OPTICAL_System')

class fold_name(path_name):
    '''
    '''
    LOG_ROOT_FOLDER = os.path.join(path_name.CALIBRATION_ROOT_FOLDER, 'LOGData/mylog')
    IFFUNCTIONS_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER, 'IFFunctions')
