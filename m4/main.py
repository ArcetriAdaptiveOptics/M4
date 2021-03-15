'''
Authors
  - C. Selmi: written in 2020

List of contents:

Functions for tower alignment
+++++++++++++++++++++++++++++
- :func:`ott_alignment_calibration`
- :func:`ott_alignment`
- :func:`m4_alignment_calibration`
- :func:`m4_alignment`
- :func:`rotation_and_optical_axis_alignment`

Functions for noise measurements
++++++++++++++++++++++++++++++++
- :func:`opto_mech_disturbances_acquisition`
- :func:`stability_vibrations`
- :func:`spectrumFromData`
- :func:`convection_moise`
- :func:`piston_noise`

PT sensors
++++++++++
- :func:`PT_calibration`
- :func:`analyzer_PT_meas`

'''

import os
import time
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.configuration import config
from m4.noise_functions import Noise
from m4.alignment import Alignment
from m4.ground import logger_set_up as lsu
from m4.configuration import start
from m4.configuration.ott_parameters import OttParameters, OpcUaParameters


def start_log(logging_level):
    """
    Parameters
    ----------
    logging_level: int
                    Warning = 30, Info = 20, Debug = 10, Notset = 0

    """
    file_path = config.fold_name.LOG_ROOT_FOLDER
    lsu.set_up_logger(file_path, logging_level)
    return file_path


####### Allineamento Torre ########

ott, interf = start.create_ott()
a = Alignment(ott)

def ott_alignment_calibration(n_frames, commandAmpVector, nPushPull, old_or_new, move):
    '''
    Parameters
    ----------------
            command_amp_vector: numpy array
                                  vector containing the movement values
                                  of the 5 degrees of freedom
            n_push_pull: int
                        number of push pull for each degree of freedom
            move: int
                1 to move the tower
                other to show matrix delta command
            old_or_new: int
                        0 for new (mixed), 1 for old (not mixed)

    Returns
    -------
            tt_tower: string
                    calibration measurement
    '''
    print('PAR + RM calibration')
    if move == 1:
        tt_tower = a.ott_calibration(n_frames, commandAmpVector, nPushPull, old_or_new, 0)
    #mask_index = 3 per il simulatore  e 0 per la mott
        return tt_tower
    else:
        mat, cmdList = a._cal._createCommandMatrix(0, commandAmpVector, old_or_new)
        plt.clf()
        plt.imshow(mat, origin='lower')
        plt.colorbar()
        return mat

def ott_alignment(tt_tower, n_images, move=1, intMatModesVector=None, commandId=None):
    '''
    Parameters
    ----------
    tt_tower: string
            calibration measurement to use for alignment
    n_images: int
            number of interferometers frames
    move: int
        1 to move the tower
        other to show commands
    Other Parameters
    ----------
    intMatModesVecor: numpy array
                    None is equal to np.array([0,1,2,3,4,5])
                    for tip, tilt, fuoco, coma, coma
    commandId: numpy array
            array containing the number of degrees of freedom to be commanded
    '''
    print('Ott alignemnt')
    par_cmd, rm_cmd = a.ott_alignment(n_images, move, intMatModesVector, commandId, tt_tower)
    print('comandi separati')
    print(par_cmd)
    print(rm_cmd)
    #check
#    if move == 1:
#	    for i in range(OttParameters.PARABOLA_DOF.size):
#	        if par_cmd[OttParameters.PARABOLA_DOF[i]] < OttParameters.parab_max_displacement[OttParameters.PARABOLA_DOF[i]]:
#	            print('ok')
#	        else:
#	            raise OSError('Par command to large')
#	    for i in range(OttParameters.RM_DOF.size):
#	        if rm_cmd[OttParameters.RM_DOF[i]] < OttParameters.rm_max_displacement[OttParameters.RM_DOF[i]]:
#	            print('ok')
#	        else:
#	            raise OSError('Rm command to large')


def m4_alignment_calibration(nFrames, commandAmpVector_ForM4Calibration=None,
                     nPushPull_ForM4Calibration=None):
    """
    Other Parameters
    ----------------
            commandAmpVector_ForM4Calibration: numpy array
                                            amplitude to be applied to m4
            nPushPull_ForM4Calibration: int
                                        number of push pull for m4 dof
            nFrames = int
                    frames for 4D

    Returns
    -------
            tt_m4: string
                    calibration measurement
    """
    print('M4 calibration')
    if commandAmpVector_ForM4Calibration is None:
        commandAmpVector_ForM4Calibration = np.array([5.0e-06, 5.0e-06])
    if nPushPull_ForM4Calibration is None:
        nPushPull_ForM4Calibration = 3
    tt_m4, zCoefComa, comaSurface = a.m4_calibration(commandAmpVector_ForM4Calibration,
                                                     nPushPull_ForM4Calibration, 5, nFrames)
    return tt_m4

def m4_alignment(tt_m4):
    '''
    Parameters
    ----------
    tt_m4: string
            calibration measurement to use for alignment
    '''
    print('M4 alignment')
    zCoefComa = a._readZcoef(tt_m4)
    cmd_m4 = a.m4_alignment(zCoefComa, tt_m4)
    print(cmd_m4)
    #check
    #applicare comando
    for i in range(OttParameters.M4_DOF.size):
        if cmd_m4[OttParameters.M4_DOF[i]] < OttParameters.m4_max_displacement[OttParameters.M4_DOF[i]]:
            print('ok')
        else:
            raise OSError('Command to large')
    #a._write_m4(cmd_m4)
    return cmd_m4

def rotation_and_optical_axis_alignment(start_point, end_point, n_points):
    '''
    Parameters
    ----------
            start_point: int
                        value of start angle
            end_point: int
                        value of end angle
            n_points:int
                    number of images desired

    Returns
    -------
        ro: object
            rotation_and_optical_axis_alignment class object
        tt: strig
            tracking number of measurement
    '''
    from m4.utils.rotation_and_optical_axis_alignment import RotOptAlign
    ro = RotOptAlign(ott)

    tt = ro.image_acquisition(start_point, end_point, n_points)

    centro, axs, raggio = ro.data_analyzer(tt)
    print(centro, axs, raggio)
    #le immagini le fa l'analyzer
    return ro, tt



######### Misure di noise ##########
def _path_noise_results(data_file_path):
    results_path = os.path.join(config.path_name.OUT_FOLDER, 'Noise')
    x = data_file_path.split("/")
    dove = os.path.join(results_path, x[len(x)-2])
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
            temp = np.tile(vec, np.int(i/2))
        elif i %2 == 1:
            #dispari
            k = i-2
            if k == 1:
                temp_pari = vec
            else:
                temp_pari = np.tile(vec, np.int((i-1)/2))
            temp = np.append(temp_pari, 1)
        template_list.append(temp)
    return template_list

def stability_vibrations(data_file_path, numbers_array, tidy_or_shuffle):
    '''
    Parameters
    ----------
        data_file_path: string
                        measurement data folder
        numbers_array: numpy array
                    vector containing integers numbers for
                    template creation
        tidy_or_shuffle: int
                        0 for tidy, 1 for shuffle
    '''
    print('Noise analysis using template')
    n = Noise()
    dove = _path_noise_results(data_file_path)
    template_list = _createTemplateList(numbers_array)

    tt_list = []
    for temp in template_list:
        tt = n.noise_analysis_from_hdf5_folder(data_file_path, tidy_or_shuffle,
                                               temp)
        time.sleep(1)
        tt_list.append(tt)

    fits_file_name = os.path.join(dove, 'trackingnumbers_%d.txt' %tidy_or_shuffle)
    file = open(fits_file_name, 'w+')
    file.write('Tidy or shuffle = %d \n' %tidy_or_shuffle)
    for tt in tt_list: 
        file.write('%s \n' %tt)
    file.close()

    rms_medio, quad_medio, n_temp = n.different_template_analyzer(tt_list)
    pyfits.writeto(os.path.join(dove, 'rms_vector_%d.fits' %tidy_or_shuffle), rms_medio, overwrite=True)
    pyfits.writeto(os.path.join(dove, 'tiptilt_vector_%d.fits' %tidy_or_shuffle), quad_medio, overwrite=True)
    pyfits.writeto(os.path.join(dove, 'n_temp_vector_%d.fits' %tidy_or_shuffle), n_temp, overwrite=True)

    plt.clf()
    plt.plot(n_temp, rms_medio, '-o', label= 'rms_medio'); plt.xlabel('n_temp')
    plt.legend()
    name = os.path.join(dove, 'rms_ntemp_%d.png' %tidy_or_shuffle)
    if os.path.isfile(name):
        os.remove(name)
    plt.savefig(name)
    plt.figure()
    plt.plot(n_temp, quad_medio, '-o', label= 'tip_tilt'); plt.xlabel('n_temp')
    plt.legend()
    name = os.path.join(dove, 'tiptilt_ntemp_%d.png' %tidy_or_shuffle)
    if os.path.isfile(name):
        os.remove(name)
    plt.savefig(name)
#     plt.figure()
#     plt.plot(freq, np.absolute(spe), '-o'); plt.xlabel('Freq[HZ]');
#     plt.ylabel('|FFT(sig)|'); plt.title('tip_tilt_%d' %tidy_or_shuffle)
#     plt.savefig(os.path.join(dove, 'tiptilt_spectrum_%d.png' %tidy_or_shuffle))
    return

def spectrumFromData(data_file_path):
    '''
    Parameters
    ----------
        data_file_path: string
                        measurement data folder
    '''
    print('Spectrum analysis')
    n = Noise()
    dove = _path_noise_results(data_file_path)

    tip, tilt = n._spectrumAllData(data_file_path)
    spe_tip, freq_tip = n._fft(tip)
    spe_tilt, freq_tilt = n._fft(tilt)

    plt.clf()
    plt.plot(freq_tip, np.absolute(spe_tip), 'o'); plt.xlabel('Freq[HZ]')
    plt.ylabel('|FFT(sig)|'); plt.title('tip_spectrum')
    name = os.path.join(dove, 'tip_spectrum.png')
    if os.path.isfile(name):
        os.remove(name)
    plt.savefig(name)
    plt.figure()
    plt.plot(freq_tilt, np.absolute(spe_tilt), 'o'); plt.xlabel('Freq[HZ]')
    plt.ylabel('|FFT(sig)|'); plt.title('tilt_spectrum')
    name = os.path.join(dove, 'tilt_spectrum.png')
    if os.path.isfile(name):
        os.remove(name)
    plt.savefig(name)

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
#     tau_c = 30 * (1/27.58)
#     epsilon_c = 2 * tau_c * np.sqrt(n_meas)
#     fits_file_name = os.path.join(dove, 'epsilon_c.txt')
#     file = open(fits_file_name, 'w+')
#     file.write('Epsilon_c = %e' %epsilon_c)
#     file.close()
    return

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
    '''
    Parameters
    ----------
        n_meas: int
            number of measurement to store

    Returns
    -------
        dove: string
            data file path of measurement
    '''
    from m4.ground import tracking_number_folder
    from opcua import Client
    server = OpcUaParameters.server
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
    return dove

def analyzer_PT_meas(tt):
    '''
    Parameters
    ----------
    tt: string
        tracking number folder
    '''
    #tt = '20200911_142702'
    from m4.ground import smooth_function

    folder = config.fold_name.PT_ROOT_FOLDER
    name = os.path.join(folder, tt)
    list = os.listdir(name)
    list.sort()

    matrix = np.zeros((len(list), OpcUaParameters.num_PT_sensor))
    matrix_s = np.zeros((len(list), OpcUaParameters.num_PT_sensor))

    i = 0
    for t in list:
        hduList = pyfits.open(os.path.join(name, t))
        temp = hduList[0].data
        matrix[i,:] = temp
        i = i+1

    matrixDiff = matrix - matrix[0,:]
    i=0
    for i in range(OpcUaParameters.num_PT_sensor):
        ss = smooth_function.smooth(matrixDiff[:,i],9)
        matrix_s[:,i] = ss
        i=i+1

    t = np.arange(0,2*len(list),2)
    plt.plot(t, matrix_s/100)
    plt.xlabel('Time [s]'); plt.ylabel('Temperature [C]'); 
    plt.title('PT Calibration')
