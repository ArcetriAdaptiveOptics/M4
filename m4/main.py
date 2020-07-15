'''
@author: cs
'''
import os
import numpy as np
from matplotlib import pyplot as plt
from m4.ground.configuration import Configuration
from m4.noise_functions import Noise
from m4.alignment import Alignment
from m4.ground import logger_set_up as lsu
from m4.configuration import start

def start_log(logging_level):
    file_path = Configuration.LOG_ROOT_FOLDER
    lsu.set_up_logger(file_path, logging_level)
#     file_path = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, 'LogIFF')
#     lsu.set_up_logger(file_path, logging_level)

ott = start.create_ott()
a = Alignment(ott)
def ott_alignement_calibration(commandAmpVector=None, nPushPull=None):
    '''
    Parameters
    ----------
            command_amp_vector: numpy array
                                  vector containing the movement values
                                  of the 5 degrees of freedom
            n_push_pull: int
                        number of push pull for each degree of freedom
            commandAmpVector_ForM4Calibration: numpy array
                                            amplitude to be applied to m4
            nPushPull_ForM4Calibration: int
                                        number of push pull for m4 dof

    Returns
    -------
            par_cmd: numpy array
                    vector of command to apply to PAR dof
            rm_cmd: numpy array
                    vector of command to apply to RM dof
            m4_cmd: numpy array
                    vector of command to apply to M4 dof
    '''
    a._moveSegmentView(0.75, 90.)
    a._moveRM(0.6)
    tt_tower = a.ott_calibration(commandAmpVector, nPushPull, 3)
    return tt_tower

def ott_alignement(tt_tower):
    par_cmd, rm_cmd = a.ott_alignement(tt_tower)
    #check
    #applicare comando (separare l'aplycmd e decidere dove metterlo)


def m4_alignement_calibration(commandAmpVector_ForM4Calibration,
                     nPushPull_ForM4Calibration):
    a._moveSegmentView(0.75, 90.)
    a._moveRM(0.6)
    tt_m4, zCoefComa, comaSurface = a.m4_calibration(commandAmpVector_ForM4Calibration,
                                                     nPushPull_ForM4Calibration, 5)
    return tt_m4

def m4_alignement(tt_m4):
    #rileggere valore coma
    cmd_m4 = a.m4_alignement(zCoefComa, tt_m4)
    #check
    #applicare comando
    return cmd_m4

###

def opto_mech_disturbances():
    #acquisizione di immagini con un certo criterio
    #definire una serie di intervalli per prendere misure
    data_file_path = os.path.join(Noise._storageFolder(), 'hdf5')
    return data_file_path

def stability_vibrations(data_file_path, template_list, tidy_or_shuffle):
    '''
    Parameters
    ----------
        data_file_path: string
                        measurement data folder
        template_list: list
                     list of template to use in the analysis
        tidy_or_shuffle: int
                        0 for tidy, 1 for shuffle

    Returns
    -------
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

def convection_noise(data_file_path, tau_vector):
    '''
    Parameters
    ----------
        data_file_path: string
                        measurement data folder
        tau_vector: numpy array
                    vector of tau to use
    '''
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
    '''
    Parameters
    ----------
        data_file_path: string
                        measurement data folder
    '''
    n = Noise()
    mean, time = n.piston_noise(data_file_path)
    plt.plot(time, mean); plt.xlabel('time[s]'); plt.ylabel('mean_image')

#PROCEDURE OTT#
def RS_verification():
    #caliball

# Mount the reference sphere in front on the RS
# Align the reference sphere to the RS
# Acquire a measurement
# Rotate the reference sphere
# Acquire a measurement
# Repeat points 4-5 until the residual measurement noise is below Test Pass Criteria
# Average the result to obtain the RS cavity
    pass

def PAR_verification():
# Mount the PAR
# Mount the INT
# Mount the CGH in front of the interferometer
# Align the CGH to the INT
# Align the INT+CGH to the PAR
# Acquire the measurement
# Post process the phasemaps to average WFE and WFE noise
    pass

def LAI_verification():
# To null the fringes:
#
# Mount the PAR
# Mount the RM and place it at the PAR center
# Remove the reference sphere
# Place a pinhole at the L1 focus
# Align the RM to have beam and return passing through the pinhole
# Align the RM to null the tilt fringes on the interferometer screen
# Align the PM and RM together to null coma on the interferometer
#
#
# To center the pupil:
#
# Identify the markers of the PM and its center.
# Adjust the pupil relay table to move it to the CCD center. Move table XY.
# Adjust the table tip/tilt to null the fringes.
    pass

def SAI_verification():
# Mount the flat reference mirror in front on the SAI. Adjust the flat tip/tilt
# Align the flat to the SAI
# Acquire a measurement
# Repeat point 3 until the residual measurement noise is below the Test Pass Criteria
# Average the result to obtain the SAI cavity WFE
    pass

def SPL_verification():
    #test 0087
    pass

# prendere ispirazione dal sito adoptica
def acquire_IFFunction():
    pass
