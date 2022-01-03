import numpy as np
from m4.utils.SPL_data_acquirer import SplAcquirer
from m4.analyzers.SPL_data_analyzer import SplAnalyzer

def define_camera():
    import pysilico
    IPServer = 'chissà'
    port = 0000
    cam = pysilico.camera(IPServer, port)
    pass

def define_filter():
    pass

def SPL_measurement_and_analysis(camera, filter):
    meas = SplAcquirer(camera, filter)
    lambda_vector = np.arange(530,730,10)
    tt = meas.acquire(lambda_vector, exptime=0.7, mask=None)
    print(tt)
    an = SplAnalyzer()
    piston, piston_smooth = an.analyzer(tt)
    print(piston, piston_smooth)
    return tt
