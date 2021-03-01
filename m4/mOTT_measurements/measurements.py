#!/User/rm/anaconda3/bin/python3.7
#!/usr/bin/python
# ./library/mycode.py
#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
'''
Authors
  - C. Selmi: written in 2020
'''

from m4.configuration import start
import time
import os
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.configuration.ott_parameters import OpcUaParameters
from m4.ground.opc_ua_controller import OpcUaController
from m4.configuration.config import fold_name
from m4.alignment import Alignment
from m4.ground import tracking_number_folder
from m4.ground import zernike
from m4.ground.timestamp import Timestamp
from m4.utils import accelerometer_data as acc


ott, interf = start.create_ott()
a = Alignment(ott)
opc = OpcUaController()


def testCalib(commandAmpVector, n_repetition=15):
    '''
    Parameters
    ----------
    commandAmpVector: numpy array
        vector containing the movement values
        of the 5 degrees of freedom
    n_repetition: int
        number of calibration acquisition and analysis to repeat

    Returns
    ------
    ttAmpVector: numpy array
        vector of tracking numbers strings of measurements
    '''
    nPushPull = 4
    n_frames = 5
    old_or_new = 1 # not mixed
    tt_list = []
    file_name = os.path.join(fold_name.CALIBRATION_ROOT_FOLDER,
                             'TtSleepCalib.txt')
    file = open(file_name, 'w+')
    for i in range(n_repetition):
        tt_tower = a.ott_calibration(n_frames, commandAmpVector,
                                     nPushPull, old_or_new, 0)
        tt_list.append(tt_tower)
        print(tt_tower)
        file.write('%s' % tt_tower)
    file.close()

    tt_list.append(commandAmpVector)
    ttAmpVector = np.array(tt_list)
    return ttAmpVector


def opticalMonitoring(n_images, delay):
    '''
    Parameters
    ----------
    n_images: int
        number of images to acquire

    Returns
    ------
    tt: string
        tracking number of measurements
    '''
    opc = OpcUaController()
    store_in_folder = fold_name.OPD_SERIES_ROOT_FOLDER
    save = tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()

    zer_list = []
    temp_list = []
    t0 = time.time()
    for i in range(n_images):
        ti = time.time()
        dt = ti - t0
        masked_ima = interf.acq4d(1, ott)
        temp_vect = opc.get_temperature_vector()
        name = Timestamp.now() + '.fits'
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, masked_ima.data)
        pyfits.append(fits_file_name, masked_ima.mask.astype(int))

        coef, mat = zernike.zernikeFit(masked_ima, np.arange(10) + 1)
        vect = np.append(dt, coef)
        zer_list.append(vect)
        temp_list.append(temp_vect)

        fits_file_name = os.path.join(dove, 'zernike.fits')
        pyfits.writeto(fits_file_name, np.array(zer_list), overwrite=True)
        fits_file_name = os.path.join(dove, 'temperature.fits')
        pyfits.writeto(fits_file_name, np.array(temp_list), overwrite=True)

        time.sleep(delay)
    return tt


def actsRepeatability(n_meas, piston_value, n_frames):
    '''
    Parameters
    ----------
    n_meas: int
        number of measurement for the test
    piston_value: float
        relative value for the parabola piston
    n_frames: int
        number of frames for the acquisition measurement

    Returns
    ------
    tt: string
        tracking number of measurements
    '''
    store_in_folder = fold_name.REPEATABILITY_ROOT_FOLDER
    save = tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()

    piston = np.array([0, 0, piston_value, 0, 0, 0])
    pos_par = ott.parab()
    pos_rm = ott.refflat()

    par_list = []
    rm_list = []
    cube = None
    for i in range(n_meas):
        ott.parab(pos_par)
        ott.refflat(pos_rm)
        par0 = _readActs(OpcUaParameters.PAR1, OpcUaParameters.PAR2,
                         OpcUaParameters.PAR3)
        rm0 = _readActs(OpcUaParameters.RM1, OpcUaParameters.RM2,
                        OpcUaParameters.RM3)
        masked_ima0 = interf.acq4d(n_frames, ott)

        ott.parab(pos_par + piston)
        # ott.refflat(pos_rm + piston)

        par1 = _readActs(OpcUaParameters.PAR1, OpcUaParameters.PAR2,
                         OpcUaParameters.PAR3)
        rm1 = _readActs(OpcUaParameters.RM1, OpcUaParameters.RM2,
                        OpcUaParameters.RM3)
        masked_ima1 = interf.acq4d(n_frames, ott)

        ott.parab(pos_par - piston)
        # ott.refflat(pos_rm - piston)

        par2 = _readActs(OpcUaParameters.PAR1, OpcUaParameters.PAR2,
                         OpcUaParameters.PAR3)
        rm2 = _readActs(OpcUaParameters.RM1, OpcUaParameters.RM2,
                        OpcUaParameters.RM3)
        masked_ima2 = interf.acq4d(n_frames, ott)

        par = np.array([par0, par1, par2])
        rm = np.array([rm0, rm1, rm2])
        cubetto = np.ma.dstack((masked_ima0, masked_ima1, masked_ima2))
        if cube is None:
            cube = cubetto
        else:
            cube = np.ma.dstack((cube, cubetto))

        par_list.append(par)
        rm_list.append(rm)
        pyfits.writeto(os.path.join(dove, 'par.fits'),
                       np.array(par_list), overwrite=True)
        pyfits.writeto(os.path.join(dove, 'rm.fits'),
                       np.array(rm_list), overwrite=True)
        pyfits.writeto(os.path.join(dove, 'images.fits'),
                       cube.data, overwrite=True)
        pyfits.append(os.path.join(dove, 'images.fits'),
                      cube.mask.astype(int), overwrite=True)

    ott.parab(pos_par)
    ott.refflat(pos_rm)
    return tt


def scanAstigmComa(stepamp, nstep, nframes=10):  # by RB 20210117.
    '''
    Parameters
    ----------
    stepamp: int
        amplitude for tip and tilt
    nstep: int
        number of step measurements
    nframes: int
        number of interferometr's frames

    Returns
    ------
    tt: string
        tracking number of measurements
    '''
    # goal: to measure coma and astigmatism at different PAR position, spanning 500 arcsec
    store_in_folder = fold_name.CALIBRATION_ROOT_FOLDER
    save = tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()
    par2rm = -2.05
    zern_vect = []
    parpos = []
    rmpos = []
    par0 = ott.parab()
    rm0 = ott.refflat()
    n2move = np.array([3, 4])
    thedirection = np.array([-1, 1])
    n_frames_alignment = 3
    tt_for_align = '20210111_152430'

    for k in n2move:
        for v in thedirection:

            for i in range(nstep):
                par1 = par0.copy()
                parmove = stepamp * i * v
                par1[k] += parmove
                print('Moving PAR[%d] by %d' % (k, parmove))
                ott.parab(par1)
                rm1 = rm0.copy()
                rmmove = stepamp * i * v * par2rm
                rm1[k] += rmmove
                ott.refflat(rm1)
                print('Moving RM[%d] by %d' % (k, rmmove))
                par_cmd, rm_cmd = a.ott_alignment(n_frames_alignment, 1,
                                                  np.array([0, 1]),
                                                  np.array([3, 4]),
                                                  tt_for_align)
                par2 = ott.parab()
                rm2 = ott.refflat()
                masked_ima = interf.acq4d(nframes, 0)
                name = Timestamp.now() + '.fits'
                fits_file_name = os.path.join(dove, name)
                pyfits.writeto(fits_file_name, masked_ima.data)
                pyfits.append(fits_file_name, masked_ima.mask.astype(int))

                coef, mat = zernike.zernikeFit(masked_ima, np.arange(10) + 1)
                zern_vect.append(coef)
                parpos.append(par2)
                rmpos.append(rm2)

                fits_file_name = os.path.join(dove, 'zernike.fits')
                pyfits.writeto(fits_file_name, np.array(zern_vect),
                               overwrite=True)
                fits_file_name = os.path.join(dove, 'PAR_positions.fits')
                pyfits.writeto(fits_file_name, np.array(parpos), overwrite=True)
                fits_file_name = os.path.join(dove, 'RM_positions.fits')
                pyfits.writeto(fits_file_name, np.array(rmpos), overwrite=True)

    ott.parab(par0)
    ott.refflat(rm0)
    return tt


def _readActs(n1, n2, n3):
    act1 = opc.get_position(n1)
    act2 = opc.get_position(n2)
    act3 = opc.get_position(n3)
    return np.array([act1, act2, act3])


def parPistonTest(piston_value, deltapos_filepath, amp, tt_for_align):
    '''
    Parameters
    ----------
    piston_value: int
        relative value for the parabola piston
    deltapos_filepath: string
        file path for par and rm delta positions
    amp: float
        tip tilt amplitude for par and rm
    tt_for_align: string
        tracking number for alignment

    Returns
    -------
    tt: string
        tracking number of measurements
    '''
        # '/home/m4/pardeltapos.fits'
    hduList = pyfits.open(deltapos_filepath)
    deltapos = hduList[0].data
    dx = deltapos[:, 0] * amp
    dy = deltapos[:, 1] * amp
    save = tracking_number_folder.TtFolder(fold_name.PISTON_TEST_ROOT_FOLDER)
    dove, tt = save._createFolderToStoreMeasurements()
    par0 = ott.parab()
    rm0 = ott.refflat()
    n_frames_meas = 10
    n_frames_alignment = 3
    rmcoeff = -2.04

    coef_list = []
    par_list = []
    rm_list = []
    for i in range(dx.size):
        if i == 0:
            print('Iteration 0')
        else:
            print('Iteration %d' % i)
            par_new = par0.copy()
            rm_new = rm0.copy()
            par_new[3] += dx[i]
            rm_new[3] += rmcoeff * dx[i]
            par_new[4] += dy[i]
            rm_new[4] += rmcoeff * dy[i]
            ott.parab(par_new)
            ott.refflat(rm_new)
            par_cmd, rm_cmd = a.ott_alignment(n_frames_alignment, 1,
                                              np.array([0, 1]),
                                              np.array([3, 4]),
                                              tt_for_align)

        par = ott.parab()
        rm = ott.refflat()
        par_list.append(par)
        rm_list.append(rm)
        masked_ima0 = interf.acq4d(n_frames_meas, ott)
        par[2] += piston_value
        ott.parab(par)
        masked_ima1 = interf.acq4d(n_frames_meas, ott)
        par[2] -= piston_value
        ott.parab(par)
        diff = masked_ima1 - masked_ima0
        name = 'diff_%04d.fits' % i
        interf.save_phasemap(dove, name, diff)
        coef, mat = zernike.zernikeFit(diff, np.arange(10) + 1)
        coef_list.append(coef)

        fits_file_name = os.path.join(dove, 'Zernike.fits')
        pyfits.writeto(fits_file_name, np.array(coef_list), overwrite=True)
        fits_file_name = os.path.join(dove, 'PAR_Positions.fits')
        pyfits.writeto(fits_file_name, np.array(par_list), overwrite=True)
        fits_file_name = os.path.join(dove, 'RM_Positions.fits')
        pyfits.writeto(fits_file_name, np.array(rm_list), overwrite=True)
    return tt


def parTiltTest(act, val_vec):
    '''
    Parameters
    ----------
    act: int
        number of the actuator to move
    val_vec: int
        value to assignee to actuator

    Returns
    -------
    quad: numpy array
        vector of rss
    coef: numpy array
        vector of zernike coefficients of image-difference
    image0: numpy masked array
        start image
    delta_list: list
        different from start image
    '''
    image0 = interf.acq4d(10, 0)
    delta_list = []
    sum_list = []
    coef_list = []
    for i in range(val_vec.size):
        opc._setAct(act, val_vec[i])
        image = interf.acq4d(10, 0)
        delta = image - image0
        delta_list.append(delta)
        coef, mat = zernike.zernikeFit(delta, np.arange(10) + 1)
        coef_list.append(coef)
        sum = np.sqrt(coef[1] ** 2 + coef[2] ** 2)
        sum_list.append(sum)
    quad = np.array(sum_list)
    plt.plot(val_vec, quad, '-o')
    return quad, np.array(coef_list), image0, delta_list


def mappingPar(shift, n_iter, tt_for_align):
    '''
    Parameters
    ----------
    shift: numpy array
            np.array([shift_par, shift_rm, rot_angle])
    n_iter: int
            number of iterations
    tt_for_align: string
        tracking number for alignment

    Returns
    -------
    tt: string
        tracking number of measurements
    '''
    n_frames_alignment = 1
    par0 = ott.parab()
    rm0 = ott.refflat()
    delta_par = []
    delta_rm = []
    delta_object = []
    delta_object2 = []
    delta_object3 = []

    save = tracking_number_folder.TtFolder(fold_name.MAPPING_TEST_ROOT_FOLDER)
    dove, tt = save._createFolderToStoreMeasurements()
    if shift[2] != 0:
        shift[0] = shift[1] = 0
        object_to_move = 'ANGLE'
    if shift[1] != 0:
        shift[0] = shift[2] = 0
        object_to_move = 'RM'
    if shift[0] != 0:
        shift[1] = shift[0]
        shift[2] = 0
        object_to_move = 'PAR + RM'

    file = open(os.path.join(dove, 'MappingInfo.txt'), "a")
    file.write('Mapping object: ' + object_to_move + '\n')
    file.close()
    slide0 = ott.slide()
    rslide0 = ott.rslide()
    angle0 = ott.angle()
    slide = -1
    rslide = -1
    angle = -1
    for i in range(n_iter):
        if shift[0] != 0:
            slide = ott.slide(slide0 + ((i + 1) * shift[0]))
#             if slide==0:
#                 raise Exception('HOMING! PAR WIN')
        if shift[1] != 0:
            rslide = ott.rslide(rslide0 + ((i + 1) * shift[1]))
#             if rslide==0:
#                 raise Exception('HOMING! RM WIN')
        if shift[2] != 0:
            angle = ott.angle(angle0 + ((i + 1) * shift[2]))

        time.sleep(5)
        par_cmd, rm_cmd = a.ott_alignment(n_frames_alignment, 1,
                                          np.array([0, 1]), np.array([3, 4]),
                                          tt_for_align)
        image = interf.acq4d(1)
        name = 'image_%04d.fits' % i
        interf.save_phasemap(dove, name, image)
        par = ott.parab()
        rm = ott.refflat()
        delta_par.append(par - par0)
        delta_rm.append(rm - rm0)
        delta_object.append(slide - slide0)
        delta_object2.append(rslide - rslide0)
        delta_object3.append(angle - angle0)

        fits_file_name = os.path.join(dove, 'delta_slide.fits')
        pyfits.writeto(fits_file_name, np.array(delta_object), overwrite=True)
        fits_file_name = os.path.join(dove, 'delta_rslide.fits')
        pyfits.writeto(fits_file_name, np.array(delta_object2), overwrite=True)
        fits_file_name = os.path.join(dove, 'delta_PAR_positions.fits')
        pyfits.writeto(fits_file_name, np.array(delta_par), overwrite=True)
        fits_file_name = os.path.join(dove, 'delta_RM_positions.fits')
        pyfits.writeto(fits_file_name, np.array(delta_rm), overwrite=True)
        fits_file_name = os.path.join(dove, 'delta_ANGLE_positions.fits')
        pyfits.writeto(fits_file_name, np.array(delta_object3), overwrite=True)

    return tt

def alignTest(tt, n_images, perturbation_vec, pre=False):
    '''
    Parameters
    ---------
    tt: string
        tracking number for the alignment
    n_images: int
        number of frames for the interferometer
    perturbation_vec: numpy array
        np.array([par_focus, par_tilt, par_tip])

    Other Parameters
    ----------------
    pre: boolean
        if True a pre-tiptilt alignment before the 5 dof one

    Returns
    -------
    coeff_matrix: numpy array
        zernike coefficients matrix for all the images
    tt: string
        tracking number of measurements
    '''
    # Homing the system
    a.ott_alignment(n_images, 1, np.array([0,1]), np.array([3,4]), tt)
    par0 = ott.parab()
    rm0 = ott.refflat()
    # Set perturbation from perturbation_vec
    focus = perturbation_vec[0]
    tip = perturbation_vec[1]
    tilt = perturbation_vec[2]
    par1 = np.copy(par0)
    par1[2] += focus
    par1[3] += tip
    par1[4] += tilt
    rm1 = np.copy(rm0)
    rm1[3] += -tip*2.05
    rm1[4] += -tilt*2.05
    # Initialization
    zern_vec = np.arange(1, 12, 1)
    coeff = []
    # Start image, perturbation and perturbed image
    image = interf.acq4d(n_images)
    pippo = zernike.zernikeFit(image, zern_vec)
    coeff.append(pippo[0])
    ott.parab(par1)
    ott.refflat(rm1)
    image = interf.acq4d(n_images)
    pippo = zernike.zernikeFit(image, zern_vec)
    coeff.append(pippo[0])
    # TipTilt pre-alignment
    if pre is True:
        a.ott_alignment(n_images, 1, np.array([0,1]), np.array([3,4]), tt)
        image = interf.acq4d(2)
        pippo = zernike.zernikeFit(image, zern_vec)
        coeff.append(pippo[0])
    # First alignment
    a.ott_alignment(n_images, 1, np.array([0,1,2,3,4]),
                    np.array([0,1,2,3,4]), tt)
    image = interf.acq4d(2)
    pippo = zernike.zernikeFit(image, zern_vec)
    coeff.append(pippo[0])
    # Second alignment
    a.ott_alignment(n_images, 1, np.array([0,1]), np.array([3,4]), tt)
    image = interf.acq4d(n_images)
    pippo = zernike.zernikeFit(image, zern_vec)
    coeff.append(pippo[0])
    # Check image
    image = interf.acq4d(n_images)
    pippo = zernike.zernikeFit(image, zern_vec)
    coeff.append(pippo[0])

    coeff_matrix = np.array(coeff)*1e9
    parend = ott.parab()

    store_in_folder = os.path.join(fold_name.REPEATABILITY_ROOT_FOLDER,
                                   'Alignment')
    save = tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()
    fits_file_name = os.path.join(dove, 'zernike.fits')
    pyfits.writeto(fits_file_name, coeff_matrix, overwrite=True)

    print(parend-par0)
    return coeff_matrix, tt

