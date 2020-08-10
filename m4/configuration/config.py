'''
Autors
  - C. Selmi: written in 2020
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
    LOG_ROOT_FOLDER = os.path.join(path_name.CALIBRATION_ROOT_FOLDER,
                                   'LOGData/mylog')
    IFFUNCTIONS_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER,
                                           'IFFunctions')
    FLAT_ROOT_FOLD = os.path.join(path_name.OPD_DATA_FOLDER, "Flattening")
    CALIBRATION_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER,
                                           "Calibration")
    ALIGNMENT_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER, "Alignment")
    ZERNIKECOMMANDTEST_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER,
                                                  "ZernikeCommandTest")
    NOISE_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER, "Noise")
    SPL_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER, "SPL")
    CALIBALL_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER,
                                        "Caliball")
    MODESVECTOR_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER,
                                           "ModesVector")
    MODALBASE_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER,
                                         "ModalBase")
    MODALAMPLITUDE_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER,
                                              "ModalAmplitude")
    COMMANDHISTORY_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER,
                                              "CommandHistory")
    GEOTRANSFORM_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER,
                                            "GeomTransf")
    ROT_OPT_ALIGN_ROOT_FOLDER = os.path.join(path_name.OPD_DATA_FOLDER,
                                             "RotOptAlignment")
    # metterci tutti i path di tutte le classi
