'''
@author: cs
'''
import numpy as np
from matplotlib import pyplot as plt
from m4.ground.configuration import Configuration
from m4.noise_functions import Noise
from m4.alignment import Alignment
from m4.ground import logger_set_up as lsu

def start_log(logging_level):
    file_path = Configuration.LOG_ROOT_FOLDER
    lsu.set_up_logger(file_path, logging_level)
#     file_path = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, 'LogIFF')
#     lsu.set_up_logger(file_path, logging_level)

def total_alignement(commandAmpVector, nPushPull, maskIndex,
               commandAmpVector_ForM4Calibration,
               nPushPull_ForM4Calibration, maskIndex_ForM4Alignement):
    a = Alignment()
    tt_tower = a.ott_calibration(commandAmpVector, nPushPull, maskIndex)
    cmd_tower = a.ott_alignement(tt_tower)
    tt_m4, zCoefComa, comaSurface = a.m4_calibration(commandAmpVector,
                                                     nPushPull, maskIndex,
                                                     commandAmpVector_ForM4Calibration,
                                                     nPushPull_ForM4Calibration,
                                                     maskIndex_ForM4Alignement)
    cmd_m4 = a.m4_alignement(zCoefComa, tt_m4)
    return cmd_tower, cmd_m4

def opto_mech_disturbances():
    #acquisizione di immagini con un certo criterio
    #definire una serie di intervalli per prendere misure
    n = Noise()
    data_file_path = n.FrameAcquisition()
    return data_file_path

def stability_vibrations(data_file_path, template_list, tidy_or_shuffle):
    '''
    args :
        template_list = list of template to use in the analysis
        tidy_or_shuffle = (int)  0 for tidy, 1 for shuffle

    retirns:
    '''
    n = Noise()

    tt_list = []
    for temp in template_list:
        tt = n.noise_analysis_from_hdf5_folder(data_file_path, tidy_or_shuffle,
                                               temp)
        tt_list.append(tt)

    rms_medio, quad_medio, n_temp = n.different_template_analyzer(tt_list)
    plt.plot(n_temp, rms_medio, '-o', label= 'rms_medio'); plt.xlabel('n_temp')
    plt.plot(n_temp, quad_medio, '-o', label= 'tip_tilt'); plt.xlabel('n_temp')
    plt.legend()
    return tt

def convection_noise(data_file_path, tau_vector):
    n = Noise()

    rms, quad, n_meas = n.analysis_whit_structure_function(data_file_path, tau_vector)
    plt.figure()
    plt.plot(tau_vector, rms, '-o', label= 'rms_medio'); plt.xlabel('tau_vector')
    plt.plot(tau_vector, quad, '-o', label= 'tip_tilt'); plt.xlabel('tau_vector')
    plt.legend()

    #stimare tc dal grafico e usare 2*tau_c = epsilon_c / np.sqrt(n) n = 4000
    tau_c = 30 * (1/27.58)
    epsilon_c = 2 * tau_c * np.sqrt(n_meas)
    return epsilon_c

def piston_noise(data_file_path):
    #togliere il tip e til dalle immagini e plottare la media
    #dovrei vedere una variazione nel tempo
    n = Noise()
    mean, time = n.piston_noise(data_file_path)
    plt.plot(time, mean); plt.xlabel('time[s]'); plt.ylabel('mean_image')

#PROCEDURE OTT#
def caliball():
    pass

# prendere ispirazione dal sito adoptica
def acquire_IFFunction():
    pass