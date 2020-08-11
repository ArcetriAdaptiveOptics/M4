'''
@author: cs
'''
import os
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
#from m4.ground.configuration import Configuration
from m4.configuration import config
from m4.noise_functions import Noise
from m4.alignment import Alignment
from m4.ground import logger_set_up as lsu
from m4.configuration import start
from m4.configuration.ott_parameters import *

def start_log(logging_level):
    file_path = config.fold_name.LOG_ROOT_FOLDER
    lsu.set_up_logger(file_path, logging_level)
#     file_path = os.path.join(Configuration.IFFUNCTIONS_ROOT_FOLDER, 'LogIFF')
#     lsu.set_up_logger(file_path, logging_level)

ott = start.create_ott()
a = Alignment(ott)
def ott_alignment_calibration(commandAmpVector=None, nPushPull=None):
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
    if commandAmpVector is None:
        commandAmpVector = np.array([5.0e-06, 5.0e-06, 5.0e-05, 5.0e-06, 5.0e-06])
    if nPushPull is None:
        nPushPull = 3
    tt_tower = a.ott_calibration(commandAmpVector, nPushPull, 3)
    return tt_tower

def ott_alignment(tt_tower):
    par_cmd, rm_cmd = a.ott_alignment(tt_tower)
    print(par_cmd, rm_cmd)
    #check
    #applicare comando (separare l'aplycmd e decidere dove metterlo)
    for i in range(OttParameters.PARABOLA_DOF.size):
        if par_cmd[OttParameters.PARABOLA_DOF[i]] < OttParameters.parab_max_displacement[OttParameters.PARABOLA_DOF[i]]:
            lala=0
        else:
            raise OSError('Par command to large')
    for i in range(OttParameters.RM_DOF.size):
        if rm_cmd[OttParameters.RM_DOF[i]] < OttParameters.rm_max_displacement[OttParameters.RM_DOF[i]]:
            lala=1
        else:
            raise OSError('Rm command to large')

    a._write_par(par_cmd)
    a._write_rm(rm_cmd)


def m4_alignment_calibration(commandAmpVector_ForM4Calibration=None,
                     nPushPull_ForM4Calibration=None):
    a._moveSegmentView(0.75, 90.)
    a._moveRM(0.6)
    if commandAmpVector_ForM4Calibration is None:
        commandAmpVector_ForM4Calibration = np.array([5.0e-06, 5.0e-06])
    if nPushPull_ForM4Calibration is None:
        nPushPull_ForM4Calibration = 3
    tt_m4, zCoefComa, comaSurface = a.m4_calibration(commandAmpVector_ForM4Calibration,
                                                     nPushPull_ForM4Calibration, 5)
    return tt_m4

def m4_alignment(tt_m4):
    zCoefComa = a._readZcoef(tt_m4)
    cmd_m4 = a.m4_alignment(zCoefComa, tt_m4)
    #check
    #applicare comando
    for i in range(OttParameters.M4_DOF.size):
        if cmd_m4[OttParameters.M4_DOF[i]] < OttParameters.m4_max_displacement[OttParameters.M4_DOF[i]]:
            lala=0
        else:
            raise OSError('Command to large')
    a._write_m4(cmd_m4)
    return cmd_m4

###

def opto_mech_disturbances_acquisition(nFrame, name=None, produce=0):
    #acquisizione di immagini con un certo criterio
    #definire una serie di intervalli per prendere misure
    #data_file_path = os.path.join(Noise._storageFolder(), 'hdf5')
    from oaautils import i4d
    from m4.ground.timestamp import Timestamp
    tt = Timestamp.now()

    interf = i4d.I4D()
    interf.connect()
    if name is None:
    	name = tt

    interf.capture(nFrame, name= name)
    if produce == 1:
    	interf.produce(name)
    interf.disconnect()

    data_file_path = os.path.join(config.fold_name.PHASECAM_ROOT_FOLDER, name)
    return data_file_path

def _path_noise_results(data_file_path):
    results_path = config.path_name.OUT_FOLDER
    x = data_file_path.split("/")
    dove = os.path.join(results_path, x[len(x)-1])
    if os.path.exists(dove):
        dove = dove
    else:
        os.makedirs(dove)
    return dove

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
    dove = _path_noise_results(data_file_path)

    tt_list = []
    for temp in template_list:
        tt = n.noise_analysis_from_hdf5_folder(data_file_path, tidy_or_shuffle,
                                               temp)
        tt_list.append(tt)

    fits_file_name = os.path.join(dove, 'trackingnumbers_%d.txt' %tidy_or_shuffle)
    file = open(fits_file_name, 'w+')
    file.write('Tidy or shuffle = %d \n' %tidy_or_shuffle)
    for tt in tt_list: 
        file.write('%s \n' %tt)
    file.close()

    rms_medio, quad_medio, n_temp = n.different_template_analyzer(tt_list)
    spe, freq = n._fft(quad_medio)
    pyfits.writeto(os.path.join(dove, 'rms_vector_%d.fits' %tidy_or_shuffle), rms_medio)
    pyfits.writeto(os.path.join(dove, 'tiptilt_vector_%d.fits' %tidy_or_shuffle), quad_medio)
    pyfits.writeto(os.path.join(dove, 'n_temp_vector_%d.fits' %tidy_or_shuffle), n_temp)

    plt.clf()
    plt.plot(n_temp, rms_medio, '-o', label= 'rms_medio'); plt.xlabel('n_temp')
    plt.legend()
    plt.savefig(os.path.join(dove, 'rms_ntemp_%d.png' %tidy_or_shuffle))
    plt.figure()
    plt.plot(n_temp, quad_medio, '-o', label= 'tip_tilt'); plt.xlabel('n_temp')
    plt.legend()
    plt.savefig(os.path.join(dove, 'tiptilt_ntemp_%d.png' %tidy_or_shuffle))
    plt.figure()
    plt.plot(freq, np.absolute(spe), '-o'); plt.xlabel('Freq[HZ]');
    plt.ylabel('|FFT(sig)|'); plt.title('tip_tilt_%d' %tidy_or_shuffle)
    plt.savefig(os.path.join(dove, 'tiptilt_spectrum_%d.png' %tidy_or_shuffle))
    return

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
    dove = _path_noise_results(data_file_path)

    rms, quad, n_meas = n.analysis_whit_structure_function(data_file_path, tau_vector)
    pyfits.writeto(os.path.join(dove, 'rms_vector_conv.fits'), rms)
    pyfits.writeto(os.path.join(dove, 'tiptilt_vector_conv.fits'), quad)
    pyfits.writeto(os.path.join(dove, 'tau_vector.fits'), tau_vector)

    plt.clf()
    plt.plot(tau_vector, rms, '-o', label= 'rms_medio'); plt.xlabel('tau_vector')
    plt.legend()
    plt.savefig(os.path.join(dove, 'rms_tau.png'))
    plt.figure()
    plt.plot(tau_vector, quad, '-o', label= 'tip_tilt'); plt.xlabel('tau_vector')
    plt.legend()
    plt.savefig(os.path.join(dove, 'tiptilt_tau.png'))

    #stimare tc dal grafico e usare 2*tau_c = epsilon_c / np.sqrt(n) n = 4000
    tau_c = 30 * (1/27.58)
    epsilon_c = 2 * tau_c * np.sqrt(n_meas)
    fits_file_name = os.path.join(dove, 'epsilon_c.txt')
    file = open(fits_file_name, 'w+')
    file.write('Epsilon_c = %e' %epsilon_c)
    file.close()
    return epsilon_c

def piston_noise(data_file_path):
    '''
    Parameters
    ----------
        data_file_path: string
                        measurement data folder
    '''
    n = Noise()
    dove = _path_noise_results(data_file_path)

    mean, time = n.piston_noise(data_file_path)
    spe, freq = n._fft(mean)
    pyfits.writeto(os.path.join(dove, 'piston_vector.fits'), mean)
    pyfits.writeto(os.path.join(dove, 'time_vector.fits'), time)

    plt.clf()
    plt.plot(time, mean); plt.xlabel('time[s]'); plt.ylabel('mean_image')
    plt.savefig(os.path.join(dove, 'piston_noise.png'))
    plt.figure()
    plt.plot(freq, np.absolute(spe), '-o'); plt.xlabel('Freq[HZ]');
    plt.ylabel('|FFT(sig)|'); plt.title('piston_power_spectrum')
    plt.savefig(os.path.join(dove, 'piston_spectrum.png'))

def rotation_and_optical_axis_alignment():
    from m4.utils.rotation_and_optical_axis_alignment import RotOptAlign
    ro = RotOptAlign(ott)

    ott.m4(np.array([0,0,0,9.9999997e-06,0,0]))
    tt = ro.acquire_image()

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
