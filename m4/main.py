'''
@author: cs
'''
import os
from m4.ground.configuration import Configuration
from m4.noise_functions import Noise
from m4.ground import logger_set_up as lsu

def start_log(logging_level):
    file_path = Configuration.LOG_ROOT_FOLDER
    lsu.set_up_logger(file_path, logging_level)
#     file_path = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, 'LogIFF')
#     lsu.set_up_logger(file_path, logging_level)

def opto_mech_disturbances():
    pass

def stability_vibrations():
    pass

def convection_noise():
    n = Noise()
    data_file_path = n.FrameAcquisition()
    tt = n.noise_analysis_from_hdf5_folder(data_file_path, tidy_or_shuffle,
                                           template, actsVector, n_push_pull)

def piston_noise():
    pass

#PROCEDURE OTT#
def acquire_IFFunction():
    pass