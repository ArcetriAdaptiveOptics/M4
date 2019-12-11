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

    IFFUNCTIONS_ROOT_FOLDER = \
        "/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions"

    V_MATRIX_FOR_SEGMENT_ROOT_811= \
        '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions/20170630_105105/modeMatrix.fits'

    V_MATRIX_FOR_SEGMENT_ROOT_892= \
        '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/IFFunctions/20170216_123645/modeMatrix.fits'

    #Interf
    LAMBA = 2
    #Ott
    #ParabolaPupilXYRadius= np.array([128, 128, 128]) #centerX, centerY, radius
    PARABOLA_PUPIL_XYRADIUS = np.array([256, 256, 256])
    PARABOLA_DOF = 3
    RM_DOF = 2
    PIXEL_SCALE = 360.5 #PIXEL/METRI
    #DM
    M4_MECHANICAL_PUPIL_XYRADIUS = np.array([458, 458, 458]) #np.array([512, 512, 512])
    M4_OPTICAL_DIAMETER = 858
    N_SEG = 6
    N_ACTS_TOT = 5352
    N_ACT_SEG = 892
    M4_DOF = 6
    REFERENCE_ANGLE_RAD = np.pi / 3
    REFERENCE_ANGLE_DEGREES = 60
    SEGMENT_DISTANCE_FROM_CENTRE = 320
    DIAMETER_IN_PIXEL_FOR_SEGMENT_IMAGES = 512
    BIG_IMAGE_DIAMETER = 1236

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
    