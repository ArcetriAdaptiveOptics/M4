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
import numpy as np
from matplotlib import pyplot as plt
from m4.configuration import config_folder_names as config
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

def showCommandMatrixBeforeCalibration(command_amp_vector):
    from m4.utils.optical_calibration import OpticalCalibration
    cal = OpticalCalibration('nulla', 'niente')
    mat, cmdList = cal.createCmatAndCmdList(command_amp_vector,
                                            np.append(OttParameters.PARABOLA_DOF,
                                                      OttParameters.RM_DOF))
    plt.clf()
    plt.imshow(mat, origin='lower')
    plt.colorbar()
    return mat

def calibrate_PARAndRM(ott, interf, n_frames, command_amp_vector, nPushPull):
    '''
    Function to be used to calibrate parabola and reference mirror dof

    Parameters
    ----------------
    ott: object
        tower
    interf: object
        interferometer
    command_amp_vector: numpy array [mm]
                    vector containing the movement values
                    of the 5 degrees of freedom
                    [par_piston, par_tip, par_tilt, rm_tip, rm_tilt]
                note: the application of par_tip corresponds to the application of rm_tip=-2.05*par_tip
                    same for par_tilt
    n_push_pull: int
                number of push pull for each degree of freedom

    Returns
    -------
            tt_tower: string
                    calibration measurement
    '''
    c_a = OttCalibAndAlign(ott, interf)
    print('PAR + RM calibration')
    tt_tower = c_a.par_and_rm_calibrator(n_frames, command_amp_vector, nPushPull)
    return tt_tower

def showCommandForParAndRmBeforeAlignement(ott, interf, tt_cal, n_images,
                                           zernike_to_be_corrected=None, dof_command_id=None):
    from m4.utils.optical_alignment import OpticalAlignment
    al = OpticalAlignment(tt_cal, ott, interf)
    print('Calculation of the alignment command for %s' %tt_cal)
    intMat, rec, cmat = al.selectModesInIntMatAndRecConstruction(zernike_to_be_corrected, dof_command_id)

    image = interf.acquire_phasemap(n_images)
    al._intMatModesVector = zernike_to_be_corrected
    total_zernike_vector, zernike_vector_selected = al.getZernikeWhitAlignerObjectOptions(image)
    print('zernike:')
    print(zernike_vector_selected)
    M = np.dot(cmat, rec)
    cmd = - np.dot(M, zernike_vector_selected)
    par_command, rm_command = al.getReorganizatedCommandForParAndRm(cmd, dof_command_id)
    print('comandi separati')
    print(par_command)
    print(rm_command)

def align_PARAndRM(ott, interf, tt_calib, n_images,
                   zernike_to_be_corrected=None, dof_command_id=None):
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
                    None is equal to np.array([0,1,2,3,4])
                    for tip, tilt, fuoco, coma, coma
    commandId: numpy array
            array containing the number of degrees of freedom to be commanded
    '''
    move = True
    print('Ott alignemnt')
    c_a = OttCalibAndAlign(ott, interf)
    par_cmd, rm_cmd, dove = c_a.par_and_rm_aligner(move, tt_calib, n_images,
                                              zernike_to_be_corrected,
                                              dof_command_id)
    tt_align = dove.split('/')[-1]
    print('comandi separati')
    print(par_cmd)
    print(rm_cmd)
    return tt_align

#### Calibrazione ed allineamneto per m4 (in cartellaBella.m4.toImplement.ott_calibrator_and_aligner ###
def calibrate_M4():
    pass

def align_M4():
    pass

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
#spostate in m4.noise


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
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'rms_image [m]', dove, 'std.png')
    # GRAFICO SLOPE
    y = np.array(slop_list)
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'rms_slope [arcsec]', dove, 'slope.png')
    # GRAFICO DIFF PISTON
    y = np.array(diff_piston_list)
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'diff_piston [m]', dove, 'diff_piston.png')
    # GRAFICO ROC
    y = np.array(roc_list)
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'roc [m]', dove, 'roc.png')
    # GRAFICO RMS 31 MM
    y = np.array(rms31)
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'rms_31mm [m]', dove, 'rms_31mm.png')
    # GRAFICO RMS 500 MM
    y = np.array(rms500)
    plotAndSaveForReqAnalysis(x, y, 'sqrt(n_frames)', 'rms_500mm [m]', dove, 'rms_500mm.png')


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


def plotAndSaveForReqAnalysis(x, y, xlabel, ylabel, dove, image_name):
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


