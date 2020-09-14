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



####### Allineamento Torre ########

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
    print('PAR + RM calibration')
    a._moveSegmentView(0.75, 90.)
    a._moveRM(0.6)
    if commandAmpVector is None:
        commandAmpVector = np.array([5.0e-06, 5.0e-06, 5.0e-05, 5.0e-06, 5.0e-06])
    if nPushPull is None:
        nPushPull = 3
    tt_tower = a.ott_calibration(commandAmpVector, nPushPull, 3)
    return tt_tower

def ott_alignment(tt_tower):
    print('Ott alignemnt')
    par_cmd, rm_cmd = a.ott_alignment(tt_tower)
    print(par_cmd)
    print(rm_cmd)
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
    print('M4 calibration')
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
    print('M4 alignment')
    zCoefComa = a._readZcoef(tt_m4)
    cmd_m4 = a.m4_alignment(zCoefComa, tt_m4)
    print(cmd_m4)
    #check
    #applicare comando
    for i in range(OttParameters.M4_DOF.size):
        if cmd_m4[OttParameters.M4_DOF[i]] < OttParameters.m4_max_displacement[OttParameters.M4_DOF[i]]:
            lala=0
        else:
            raise OSError('Command to large')
    a._write_m4(cmd_m4)
    return cmd_m4

def rotation_and_optical_axis_alignment(tt=None):
    from m4.utils.rotation_and_optical_axis_alignment import RotOptAlign
    ro = RotOptAlign(ott)

    if tt is None:
        tt = ro.acquire_image()
    else:
        tt = tt

    centro, axs, raggio = ro.analyzer(tt)
    #le immagini le fa l'analyzer
    return centro, axs, raggio

######### Misure di noise ##########

def opto_mech_disturbances_acquisition(nFrame, name=None, produce=0):
    #acquisizione di immagini con un certo criterio
    #definire una serie di intervalli per prendere misure
    #data_file_path = os.path.join(Noise._storageFolder(), 'hdf5')
    print('Frame acquisition')
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
    results_path = os.path.join(config.path_name.OUT_FOLDER, 'Noise')
    x = data_file_path.split("/")
    dove = os.path.join(results_path, x[len(x)-1])
    if os.path.exists(dove):
        dove = dove
    else:
        os.makedirs(dove)
    return dove

def _createTemplateList(numbers_array):
    '''
    Parameters
    ----------
        numbers_array: numpy array
                    vector containing integers numbers for
                    template creation
    Returns
    -------
        template_list: list
                    list of template to use
    '''
    template_list = []
    vec = np.array([1, -1])
    for i in numbers_array:
        if i % 2 == 0:
            #pari
            k = i-2
            temp = np.tile(vec, k)
        elif i %2 == 1:
            #dispari
            k = i-2
            if k == 1:
                temp_pari = vec
            else:
                temp_pari = np.tile(vec, k-1)
            temp = np.append(temp_pari, 1)
        template_list.append(temp)
    return template_list

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
    print('Noise analysis using template')
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
#     plt.figure()
#     plt.plot(freq, np.absolute(spe), '-o'); plt.xlabel('Freq[HZ]');
#     plt.ylabel('|FFT(sig)|'); plt.title('tip_tilt_%d' %tidy_or_shuffle)
#     plt.savefig(os.path.join(dove, 'tiptilt_spectrum_%d.png' %tidy_or_shuffle))
    return

def spectrumFromData(data_file_path):
    print('Spectrum analysis')
    n = Noise()
    dove = _path_noise_results(data_file_path)

    tip, tilt = n._spectrumAllData(data_file_path)
    spe_tip, freq_tip = n._fft(tip)
    spe_tilt, freq_tilt = n._fft(tilt)

    plt.clf()
    plt.plot(freq_tip, np.absolute(spe_tip), 'o'); plt.xlabel('Freq[HZ]')
    plt.ylabel('|FFT(sig)|'); plt.title('tip_spectrum')
    plt.savefig(os.path.join(dove, 'tip_spectrum.png'))
    plt.figure()
    plt.plot(freq_tilt, np.absolute(spe_tilt), 'o'); plt.xlabel('Freq[HZ]')
    plt.ylabel('|FFT(sig)|'); plt.title('tilt_spectrum')
    plt.savefig(os.path.join(dove, 'tilt_spectrum.png'))

def convection_noise(data_file_path, tau_vector):
    '''
    Parameters
    ----------
        data_file_path: string
                        measurement data folder
        tau_vector: numpy array
                    vector of tau to use
    '''
    print('Noise analysis using tau vector')
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
    plt.plot(freq, np.absolute(spe), 'o'); plt.xlabel('Freq[HZ]');
    plt.ylabel('|FFT(sig)|'); plt.title('piston_power_spectrum')
    plt.savefig(os.path.join(dove, 'piston_spectrum.png'))



######## Sensori PT #######

def PT_calibration(n_meas):
    from m4.ground import tracking_number_folder
    from opcua import Client
    import time
    server = "opc.tcp://192.168.22.100:48050"
    client = Client(url=server)
    client.connect()

    folder = config.fold_name.PT_ROOT_FOLDER
    save = tracking_number_folder.TtFolder(folder)
    dove, tt = save._createFolderToStoreMeasurements()

    for i in range(n_meas):
        time.sleep(2)
        temp = client.get_node("ns=7;s=MAIN.i_Temperature_Sensor")
        temp_list = temp.get_value()
        temp_vector = np.array(temp_list.get_value())

        fits_file_name = os.path.join(dove, 'temperature_%04.fits' %i)
        pyfits.writeto(fits_file_name, temp_vector)

        print('Misura %04d' %i)

def analyzer_PT_meas(tt):
    from m4.ground import smooth_function

    folder = config.fold_name.PT_ROOT_FOLDER
    name = os.path.join(folder, '20200911_142702')
    list = os.listdir(name)
    list.sort()

    matrix = np.zeros((len(list), 24))
    matrix_s = np.zeros((len(list), 24))

    i = 0
    for t in list:
        hduList = pyfits.open(os.path.join(name, t))
        temp = hduList[0].data
        matrix[i,:] = temp
        i = i+1

    matrixDiff = matrix - matrix[0,:]
    i=0
    for i in range(24):
        ss = smooth_function.smooth(matrixDiff[:,i],9)
        matrix_s[:,i] = ss
        i=i+1

    t = np.arange(0,2*len(list),2)
    plt.plot(t, matrix_s/100)
    plt.xlabel('Time [s]'); plt.ylabel('Temperature [C]'); 
    plt.title('PT Calibration')


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
