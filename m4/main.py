"""
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

"""
import os
import numpy as np
from matplotlib import pyplot as plt
from m4.configuration import config_folder_names as config
from m4.ott_calibrator_and_aligner import OttCalibAndAlign
from m4.ground import logger_set_up as lsu
from m4.configuration.ott_parameters import OttParameters
import playsound
from m4.configuration.ott_parameters import Sound
from m4.ground import geo

# print("sono io")
def start_log(logging_level=20):
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
    """
    Parameters
    ----------------
    command_amp_vector: numpy array [mm, arcsec, arcsec, arcsec, arcsec]
                    vector containing the movement values
                    of the 5 degrees of freedom
                    [par_piston, par_tip, par_tilt, rm_tip, rm_tilt]
                note: the application of par_tip corresponds to the application of rm_tip=-2.05*par_tip
                    same for par_tilt
    Returns
    -------
    mat: numpy array
        matrix of command
    """
    from m4.utils.optical_calibration import OpticalCalibration

    cal = OpticalCalibration("nulla", "niente")
    mat, cmdList = cal.createCmatAndCmdList(
        command_amp_vector, np.append(OttParameters.PARABOLA_DOF, OttParameters.RM_DOF)
    )
    plt.clf()
    x_old = np.arange(mat.shape[0])
    x = ["par_piston", "par_tip", "par_tilt", "rm_tip", "rm_tilt"]
    plt.imshow(mat, origin="lower")
    plt.xticks(x_old, x, rotation=0)
    plt.ylabel("Commands")
    plt.colorbar()
    return mat


def calibrate_PARAndRM(
    ott, interf, command_amp_vector, nPushPull, n_frames=None, delay=None,tnPar = None
):
    """
    Function to be used to calibrate parabola and reference mirror dof

    Parameters
    ----------------
    ott: object
        tower
    interf: object
        interferometer
    command_amp_vector: numpy array [mm, arcsec, arcsec, arcsec, arcsec]
                    vector containing the movement values
                    of the 5 degrees of freedom
                    [par_piston, par_tip, par_tilt, rm_tip, rm_tilt]
                note: the application of par_tip corresponds to the application of rm_tip=-2.05*par_tip
                    same for par_tilt
    n_push_pull: int
                number of push pull for each degree of freedom

    Other Parameters
    ----------------
    nframes: int
        number of frames
    delay: int [s]
        delay between images

    Returns
    -------
            tt_tower: string
                    tracking number of calibration measurements
    """
    command_amp_vector = np.array(command_amp_vector)
    if n_frames is None:
        n_frames = 1
    if delay is None:
        delay = 0
    c_a = OttCalibAndAlign(ott, interf)
    print("PAR + RM calibration")
    tt_tower = c_a.par_and_rm_calibrator(command_amp_vector, nPushPull, n_frames, delay, tnPar)
    return tt_tower


def showCommandForParAndRmBeforeAlignement(
    ott,
    interf,
    tt_cal,
    n_images,
    delay,
    zernike_to_be_corrected=None,
    dof_command_id=None,
):
    """
    Parameters
    ----------
    ott: object
        tower
    interf: object
        interferometer
    tt_cal: string
            calibration measurement to use for alignment
    n_images: int
            number of interferometers frames

    Other Parameters
    ----------
    zernike_to_be_corrected: numpy array
                    None is equal to np.array([0,1,2,3,4])
                    for tip, tilt, fuoco, coma, coma
    dof_command_id: numpy array
            array containing the number of degrees of freedom to be commanded
    """
    from m4.utils.optical_alignment import OpticalAlignment

    al = OpticalAlignment(tt_cal, ott, interf)
    print("Calculation of the alignment command for %s" % tt_cal)
    intMat, rec, cmat = al.selectModesInIntMatAndRecConstruction(
        zernike_to_be_corrected, dof_command_id
    )

    image = interf.acquire_phasemap(n_images, delay)
    al._intMatModesVector = zernike_to_be_corrected
    total_zernike_vector, zernike_vector_selected = (
        al.getZernikeWhitAlignerObjectOptions(image)
    )
    print("zernike:")
    print(zernike_vector_selected)
    M = np.dot(cmat, rec)
    cmd = -np.dot(M, zernike_vector_selected)
    par_command, rm_command = al.getReorganizatedCommandForParAndRm(cmd, dof_command_id)
    print("comandi separati")
    print(par_command)
    print(rm_command)


def align_PARAndRM(
    ott,
    interf,
    tt_calib,
    zernike_to_be_corrected=None,
    dof_command_id=None,
    n_frames=None,
    doit=False,
    tnPar = None,
    delay = None
):
    """
    Parameters
    ----------
    ott: object
        tower
    interf: object
        interferometer
    tt_tower: string
            calibration measurement to use for alignment

    Other Parameters
    ----------
    zernike_to_be_corrected: numpy array
                    None is equal to np.array([0,1,2,3,4])
                    for tip, tilt, fuoco, coma, coma
    dof_command_id: numpy array
            array containing the number of degrees of freedom to be commanded
    n_frames: int
            number of interferometers frames
    delay: int [s]
                delay between images
    """
    if n_frames is None:
        n_frames = 1
    if delay is None:
        delay = 0

    move = doit
    c_a = OttCalibAndAlign(ott, interf)
    par_cmd, rm_cmd, dove = c_a.par_and_rm_aligner(
        move, tt_calib, n_frames, delay, zernike_to_be_corrected, dof_command_id, tnPar
    )
    tt_align = dove.split("/")[-1]
    print("Mixed par+rm commands")
    print(par_cmd)
    print(rm_cmd)
    if Sound.PLAY is True:
        playsound.playsound(os.path.join(Sound.AUDIO_FILE_PATH, "aligncomp.mp3"))
    return tt_align


#### Calibrazione ed allineamneto per m4 (in cartellaBella.m4.toImplement.ott_calibrator_and_aligner ###
def calibrate_M4():
    """to be implemented"""
    pass


def align_M4():
    """to be implemented"""
    pass


def getOttConfigurator(ott):
    """Get the Class to manage ott configurations"""
    from m4.ground.ott_configurations import OttConfigurations

    oc = OttConfigurations(ott)
    return oc


def spiralize(ott, npos, step):
    p = geo.spiral_pos(npos, step)
    npos = int(len(p) / 2)
    p0 = ott.parabola.getPosition()
    for i in range(npos):
        #p0 = ott.parabola.getPosition()
        p1 = p0 + np.array([0, 0, 0, p[i, 0], p[i, 1], 0])
        print("New Par command:")
        print(p1)
        ott.parabola.setPosition(p1)
        #r0 = ott.referenceMirror.getPosition()
        #r1 = r0 + np.array([0, 0, 0, -2 * p[i, 0], -2 * p[i, 1], 0])
        #print("New RM command:")
        #print(r1)
        #ott.referenceMirror.setPosition(r1)


### ROTATION FOR ALIGNMENT (ASSE OTTICO E ASSE MECCANICO)
# tolto per usare direttamente la sua classe. Vecchia funzione in CartellaBella.copie.main
### MISURE DI NOISE ###
# spostate in m4.noise


### ANALISI DEI REQUISITI ###
# spostata in m4.requirements_checker
