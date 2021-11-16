'''
Authors
  - C. Selmi: written in 2020

List of contents:

Functions for tower alignment
+++++++++++++++++++++++++++++
- :func:`showCommandMatrixBeforeCalibration`
- :func:`calibrate_PARAndRM`
- :func:`showCommandForParAndRmBeforeAlignement`
- :func:`align_PARAndRM`
- :func:`calibrate_M4`
- :func:`align_M4`

'''

import numpy as np
from matplotlib import pyplot as plt
from m4.configuration import config_folder_names as config
from m4.ott_calibrator_and_aligner import OttCalibAndAlign
from m4.ground import logger_set_up as lsu
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



### ROTATION FOR ALIGNMENT (ASSE OTTICO E ASSE MECCANICO)
#tolto per usare direttamente la sua classe. Vecchia funzione in CartellaBella.copie.main
### MISURE DI NOISE ###
#spostate in m4.noise


### ANALISI DEI REQUISITI ###
#spostata in m4.requirements_checker
