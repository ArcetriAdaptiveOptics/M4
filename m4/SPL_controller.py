'''
Authors
    - C. Selmi: written in 2021
'''
import os
import numpy as np
from m4.utils.SPL_data_acquirer import SplAcquirer
from m4.analyzers.SPL_data_analyzer import SplAnalyzer
from m4.configuration import config_folder_names as fold_name

def define_camera():
    ''' Function to use to define the camera with pysilico
    '''
    import pysilico
    IPServer = 'chissà'
    port = 0000
    cam = pysilico.camera(IPServer, port)
    return cam

def define_filter():
    ''' Function to use to define the tunable filter
    '''
    from m4.devices.tuna_filter import TunaFilt
    conf_file = os.path.join(fold_name.CONFIGURATION_ROOT_FOLDER, 'lctf.yaml')
    filter = TunaFilt(conf_file)
    filter.connect()
    return filter

def SPL_measurement_and_analysis(camera, filter):
    '''Function for SPL data acquisition and analysis
    '''
    meas = SplAcquirer(camera, filter)
    lambda_vector = np.arange(530, 730, 10)
    tt = meas.acquire(lambda_vector, exptime=0.7, mask=None)
    print(tt)
    an = SplAnalyzer()
    piston, piston_smooth = an.analyzer(tt)
    print(piston, piston_smooth)
    return tt, piston
