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
import glob
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.configuration import config_folder_names as config
from m4.noise_functions import Noise
from m4.ott_calibrator_and_aligner import OttCalibAndAlign
from m4.ground import logger_set_up as lsu
from m4.utils import req_check
from m4.configuration.ott_parameters import OttParameters


def start_log(logging_level):
    """
    Parameters
    ----------
    logging_level: int
                    Warning = 30, Info = 20, Debug = 10, Notset = 0

    """
    file_path = config.LOG_ROOT_FOLDER
    lsu.set_up_logger(file_path, logging_level)
    return file_path

####### Allineamento Torre ########


# ott, interf = start.create_ott('/mnt/data/M4/Data/SYSCONFData/Config.yaml')
# a = Alignment(ott, interf)


def ott_calibrator(ott, interf, n_frames, command_amp_vector, nPushPull, move):
    '''
    Parameters
    ----------------
    ott: object
        tower
    interf: object
        interferometer
    command_amp_vector: numpy array
                          vector containing the movement values
                          of the 5 degrees of freedom
    n_push_pull: int
                number of push pull for each degree of freedom
    move: boolean
        True to move the tower
        other to show matrix delta command


    Returns
    -------
            tt_tower: string
                    calibration measurement
    '''
    c_a = OttCalibAndAlign(ott, interf)
    print('PAR + RM calibration')
    if move is True:
        tt_tower = c_a.ott_calibration(n_frames, command_amp_vector, nPushPull)
        return tt_tower
    else:
        mat, cmdList = c_a._cal._createCmatAndCmdList(command_amp_vector,
                                                      np.append(OttParameters.PARABOLA_DOF,
                                                                OttParameters.RM_DOF))
        plt.clf()
        plt.imshow(mat, origin='lower')
        plt.colorbar()
        return mat


def ott_aligner(ott, interf, tt_calib, n_images, move, zernike_to_be_corrected=None, dof_command_id=None):
    '''
    Parameters
    ----------
    tt_tower: string
            calibration measurement to use for alignment
    n_images: int
            number of interferometers frames
    move: boolean
        True to move the tower
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
    c_a = OttCalibAndAlign(ott, interf)
    par_cmd, rm_cmd = c_a.ott_alignment( move, tt_calib, n_images,
                      zernike_to_be_corrected, dof_command_id)
    print('comandi separati')
    print(par_cmd)
    print(rm_cmd)

#### Calibrazione ed allineamneto per m4 (in cartellaBella.m4.toImplement.ott_calibrator_and_aligner ###


### Rotation for alignment
def rotation_and_optical_axis_alignment(ott, interf, start_point, end_point, n_points):
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
    ro = RotOptAlign(ott, interf)

    tt = ro.image_acquisition(start_point, end_point, n_points)

    centro, axs, raggio = ro.data_analyzer(tt)
    print(centro, axs, raggio)
    # le immagini le fa l'analyzer
    return ro, tt


######### Misure di noise ##########
def _path_noise_results(data_file_path, h5_or_fits=None):
    ''' Function to get tt'''
    results_path = os.path.join(config.OUT_FOLDER, 'Noise')
    x = data_file_path.split("/")
    if h5_or_fits is None:
        dove = os.path.join(results_path, x[len(x) - 2])
    else:
        dove = os.path.join(results_path, x[len(x) - 1])
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
            # pari
            k = i - 2
            temp = np.tile(vec, np.int(i / 2))
        elif i % 2 == 1:
            # dispari
            k = i - 2
            if k == 1:
                temp_pari = vec
            else:
                temp_pari = np.tile(vec, np.int((i - 1) / 2))
            temp = np.append(temp_pari, 1)
        template_list.append(temp)
    return template_list


def noise_vibrations(data_file_path, numbers_array, tidy_or_shuffle):
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

    fits_file_name = os.path.join(dove, 'trackingnumbers_%d.txt' % tidy_or_shuffle)
    file = open(fits_file_name, 'w+')
    file.write('Tidy or shuffle = %d \n' % tidy_or_shuffle)
    for tt in tt_list:
        file.write('%s \n' % tt)
    file.close()

    rms_medio, quad_medio, n_temp, ptv_medio = n.different_template_analyzer(tt_list)
    pyfits.writeto(os.path.join(dove, 'rms_vector_%d.fits' % tidy_or_shuffle), rms_medio, overwrite=True)
    pyfits.writeto(os.path.join(dove, 'tiptilt_vector_%d.fits' % tidy_or_shuffle), quad_medio, overwrite=True)
    pyfits.writeto(os.path.join(dove, 'n_temp_vector_%d.fits' % tidy_or_shuffle), n_temp, overwrite=True)
    pyfits.writeto(os.path.join(dove, 'ptv_%d.fits' % tidy_or_shuffle), ptv_medio, overwrite=True)

    tt = data_file_path.split('/')[-2]
    plt.clf()
    # WFE = 2*rms_medio
    plt.plot(n_temp, rms_medio * 1e9, '-o')
    plt.xlabel('n_temp')
    plt.ylabel('rms [nm]')
    plt.title('%s' % tt)
    plt.grid()
    name = os.path.join(dove, 'rms_ntemp_%d.png' % tidy_or_shuffle)
    if os.path.isfile(name):
        os.remove(name)
    plt.savefig(name)

    plt.figure()
    plt.plot(n_temp, quad_medio * 1e9, '-o'); plt.xlabel('n_temp')
    plt.ylabel('TipTilt [nm]')
    plt.title('%s' % tt)
    plt.grid()
    name = os.path.join(dove, 'tiptilt_ntemp_%d.png' % tidy_or_shuffle)
    if os.path.isfile(name):
        os.remove(name)
    plt.savefig(name)

    plt.figure()
    plt.plot(n_temp, ptv_medio * 1e9, '-o'); plt.xlabel('n_temp')
    plt.ylabel('PtV [nm]')
    plt.title('%s' % tt)
    plt.grid()
    name = os.path.join(dove, 'ptv_ntemp_%d.png' % tidy_or_shuffle)
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

    Other Parameters
    ----------------
        h5_or_fits: if it is none the h5 data analysis is performed
    '''
    last_name = data_file_path.split('/')[-1]
    if last_name == 'hdf5':
        h5_or_fits = None
    else:
        h5_or_fits = 7

    print('Noise analysis using tau vector')
    n = Noise()
    dove = _path_noise_results(data_file_path, h5_or_fits)

    rms, quad, n_meas = n.analysis_whit_structure_function(data_file_path,
                                                           tau_vector,
                                                           h5_or_fits)
    pyfits.writeto(os.path.join(dove, 'rms_vector_conv.fits'), rms,
                   overwrite=True)
    pyfits.writeto(os.path.join(dove, 'tiptilt_vector_conv.fits'), quad,
                   overwrite=True)
    pyfits.writeto(os.path.join(dove, 'tau_vector.fits'), tau_vector,
                   overwrite=True)

    rms_nm = rms * 1e9
    if h5_or_fits is None:
        x = tau_vector * (1 / 27.58)
        param = [5, 0.5, 32]
        try:

            pp, fit = _curvFit(param, x, rms_nm)
            decorr_time = 1 / pp[0] + pp[1]
        except:
            pp = np.array([0, 0, rms[-1] * 1e9])
            decorr_time = -1
            fit = rms_nm.copy() * 0
        plt.clf()
        plt.plot(x, rms * 1e9, '-o', label='meas')
        plt.xlabel('time [s]')
        plt.ylabel('rms [nm]')
        plt.plot(x, fit, '-', label='fit')
        plt.grid()
        plt.plot([x[0], x[-1]], [pp[2], pp[2]], '--r', linewidth=3,
                 label='%.2f [nm]' % pp[2])
#         plt.plot(decorr_time, _funFit(decorr_time,*pp), 'og',
#                  label='Dec time = %d [s]' %np.round(decorr_time))
        plt.legend()
        tt = dove.split('/')[-1]
        plt.title('%s' % tt)
        name = os.path.join(dove, 'rms_tau.png')
        if os.path.isfile(name):
            os.remove(name)
        plt.savefig(name)
        return pp[2], decorr_time
    else:
        time_diff = _time_for_plot(data_file_path)
        x = tau_vector * time_diff
        plt.clf()
        plt.plot(x, rms * 1e9, '-o', label='time_diff = %d' % time_diff)
        plt.xlabel('time [s]')
        plt.ylabel('rms [nm]')
        plt.grid()
        plt.legend()
        tt = dove.split('/')[-1]
        plt.title('%s' % tt)
        name = os.path.join(dove, 'rms_tau.png')
        if os.path.isfile(name):
            os.remove(name)
        plt.savefig(name)
    # stimare tc dal grafico e usare 2*tau_c = epsilon_c / np.sqrt(n) n = 4000
#     tau_c = 30 * (1/27.58)
#     epsilon_c = 2 * tau_c * np.sqrt(n_meas)
#     fits_file_name = os.path.join(dove, 'epsilon_c.txt')
#     file = open(fits_file_name, 'w+')
#     file.write('Epsilon_c = %e' %epsilon_c)
#     file.close()


def _time_for_plot(stab_path):
    listtot = glob.glob(os.path.join(stab_path, '*.fits'))
    listtot.sort()
    aa = listtot[0].split('/')
    t0 = aa[-1].split('_')[1].split('.')[0]
    bb = listtot[1].split('/')
    t1 = bb[-1].split('_')[1].split('.')[0]

    hs = float(t0[0: 2]) * 3600
    ms = float(t0[2: 4]) * 60
    s = float(t0[4::])
    t0s = hs + ms + s
    hs = float(t1[0: 2]) * 3600
    ms = float(t1[2: 4]) * 60
    s = float(t1[4::])
    t1s = hs + ms + s

    time_diff = t1s - t0s
    return time_diff


def _funFit(x, a, b, c):
    fun = -np.exp(-a * (x - b)) + c
    return fun


def _curvFit(param, x, rms_nm):
    from scipy.optimize import curve_fit
    pp, pcov = curve_fit(_funFit, x, rms_nm, param)
    fit = _funFit(x, *pp)
    return pp, fit


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


### ANALISI DEI REQUISITI ###
def analysis_req(data_file_path, zernike_vector_to_subtract, step=None, offset=None):
    ''' Simultaneous analysis of noise requirements for a tn

    Parameters
    ----------
    path: string
        total path for data analysis

    Other Parameters
    ----------------
    offset: if it is None data analysis is made by split n_images in two
    '''
    last_name = data_file_path.split('/')[-1]
    if last_name == 'hdf5':
        tt = data_file_path.split('/')[-2]
    else:
        tt = data_file_path.split('/')[-1]

    results_path = os.path.join(config.OUT_FOLDER, 'Req')
    dove = os.path.join(results_path, tt)
    if os.path.exists(dove):
        dove = dove
    else:
        os.makedirs(dove)
    fits_file_name = os.path.join(dove, 'info.txt')
    file = open(fits_file_name, 'w+')
    if offset is None:
        file.write('Data produced without offset optic image')
    else:
        file.write('Data produced with offset optic image')
    file.close()

    print('Creating cube 50')
    image50 = req_check.robustImageFromDataSet(50, data_file_path, zernike_vector_to_subtract, offset)
    print('Creating cube 100')
    image100 = req_check.robustImageFromDataSet(100, data_file_path, zernike_vector_to_subtract, offset)
    print('Creating cube 300')
    image300 = req_check.robustImageFromDataSet(300, data_file_path, zernike_vector_to_subtract, offset)
#     print('Creating cube 600')
#     image600 = req_check.robustImageFromDataSet(600, data_file_path, offset)

    image_list = [image50, image100, image300]  # , image600]
    slop_list, diff_piston_list, roc_list, rms31, rms500 = fromImagesToReq(image_list, None, step)

    x = np.array([50, 100, 300])  # ,600])
    # GRAFICO STD IMAGES
    y = np.array([image50.std(), image100.std(), image300.std()])  # ,image600.std()])
    plotAndSave(x, y, 'sqrt(n_frames)', 'rms_image [m]', dove, 'std.png')
    # GRAFICO SLOPE
    y = np.array(slop_list)
    plotAndSave(x, y, 'sqrt(n_frames)', 'rms_slope [arcsec]', dove, 'slope.png')
    # GRAFICO DIFF PISTON
    y = np.array(diff_piston_list)
    plotAndSave(x, y, 'sqrt(n_frames)', 'diff_piston [m]', dove, 'diff_piston.png')
    # GRAFICO ROC
    y = np.array(roc_list)
    plotAndSave(x, y, 'sqrt(n_frames)', 'roc [m]', dove, 'roc.png')
    # GRAFICO RMS 31 MM
    y = np.array(rms31)
    plotAndSave(x, y, 'sqrt(n_frames)', 'rms_31mm [m]', dove, 'rms_31mm.png')
    # GRAFICO RMS 500 MM
    y = np.array(rms500)
    plotAndSave(x, y, 'sqrt(n_frames)', 'rms_500mm [m]', dove, 'rms_500mm.png')


def fromImagesToReq(image_list, pscale=None, step=None, n_patches=None):
    slop_list = []
    diff_piston_list = []
    roc_list = []
    rms31 = []
    rms500 = []
    print(pscale)
    for image in image_list:
        print('Producing slope')
        slop_list.append(req_check.test242(image, pscale))
        print('Producing differential piston')
        diff_piston_list.append(req_check.diffPiston(image))
        print('Producing roc')
        roc_list.append(req_check.test283(image, pscale, step))
        print('Producing rms31')
        rms31.append(req_check.test243(image, 0.015, pscale, step, n_patches))
        print('Producing rms51')
        rms500.append(req_check.test243(image, 0.1, pscale, step, n_patches))
    return slop_list, diff_piston_list, roc_list, rms31, rms500


def plotAndSave(x, y, xlabel, ylabel, dove, image_name):
    plt.figure(figsize=(10, 6))
    plt.plot(np.sqrt(x), y, '-o')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    tt = dove.split('/')[-1]
    plt.title('%s' % tt)
    name = os.path.join(dove, image_name)
    if os.path.isfile(name):
        os.remove(name)
    plt.savefig(name)


