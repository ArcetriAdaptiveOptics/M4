'''
@author: cs
'''
import numpy as np


class Configuration():
    '''
    Constants of system configuration
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
    M4_PUPIL_XYRADIUS = np.array([512, 512, 512])
    N_SEG = 6
    N_ACTS_TOT = 5352
    N_ACT_SEG = 892
    M4_DOF = 6
    REFERENCE_ANGLE = 60
    SEGMENT_DISTANCE_FROM_CENTRE = 28

    #per l'applicazione
    def set_up_logger(self, file_path, logging_level):
        """ Set up logger for the application """
        import logging
        import logging.handlers
        FORMAT = '%(asctime)s %(levelname)s %(name)s %(message)s'
        f = logging.Formatter(fmt = FORMAT)
        handler = logging.handlers.RotatingFileHandler(file_path,
                                                       encoding='utf8',
                                                       maxBytes=10000,
                                                       backupCount=3)
        root_logger = logging.getLogger()
        root_logger.setLevel(logging_level)
        handler.setFormatter(f)
        handler.setLevel(logging_level)
        root_logger.addHandler(handler)
        handler.doRollover()
    