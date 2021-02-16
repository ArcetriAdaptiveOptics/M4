#!/User/rm/anaconda3/bin/python3.7
#!/usr/bin/python
# ./library/mycode.py
#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
'''
Authors
  - C. Selmi: written in 2020
'''

from m4.configuration import start as _start
import time as _time
import os as _os
import numpy as _np
from astropy.io import fits as _pyfits
from matplotlib import pyplot as _plt
from m4.configuration.ott_parameters import OpcUaParameters as _OpcUaParameters
from m4.ground.opc_ua_controller import OpcUaController as _OpcUaController
from m4.configuration.config import fold_name as _fold_name
from m4.alignment import Alignment as _Alignment
from m4.ground import tracking_number_folder as _tracking_number_folder
from m4.ground import zernike as _zernike
from m4.ground.timestamp import Timestamp as _Timestamp


_ott, _interf = _start.create_ott()
_a = _Alignment(_ott)
_opc = _OpcUaController()


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
    file_name = _os.path.join(_fold_name.CALIBRATION_ROOT_FOLDER,
                             'TtSleepCalib.txt')
    file = open(file_name, 'w+')
    for i in range(n_repetition):
        tt_tower = _a.ott_calibration(n_frames, commandAmpVector,
                                     nPushPull, old_or_new, 0)
        tt_list.append(tt_tower)
        print(tt_tower)
        file.write('%s' % tt_tower)
    file.close()

    tt_list.append(commandAmpVector)
    ttAmpVector = _np.array(tt_list)
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
    opc = _OpcUaController()
    store_in_folder = _fold_name.OPD_SERIES_ROOT_FOLDER
    save = _tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()

    zer_list = []
    temp_list = []
    t0 = _time.time()
    for i in range(n_images):
        ti = _time.time()
        dt = ti - t0
        masked_ima = _interf.acq4d(1, _ott)
        temp_vect = opc.get_temperature_vector()
        name = _Timestamp.now() + '.fits'
        fits_file_name = _os.path.join(dove, name)
        _pyfits.writeto(fits_file_name, masked_ima.data)
        _pyfits.append(fits_file_name, masked_ima.mask.astype(int))

        coef, mat = _zernike.zernikeFit(masked_ima, _np.arange(10) + 1)
        vect = _np.append(dt, coef)
        zer_list.append(vect)
        temp_list.append(temp_vect)

        fits_file_name = _os.path.join(dove, 'zernike.fits')
        _pyfits.writeto(fits_file_name, _np.array(zer_list), overwrite=True)
        fits_file_name = _os.path.join(dove, 'temperature.fits')
        _pyfits.writeto(fits_file_name, _np.array(temp_list), overwrite=True)

        _time.sleep(delay)
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
    store_in_folder = _fold_name.REPEATABILITY_ROOT_FOLDER
    save = _tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()

    piston = _np.array([0, 0, piston_value, 0, 0, 0])
    pos_par = _ott.parab()
    pos_rm = _ott.refflat()

    par_list = []
    rm_list = []
    cube = None
    for i in range(n_meas):
        _ott.parab(pos_par)
        _ott.refflat(pos_rm)
        par0 = _readActs(_OpcUaParameters.PAR1, _OpcUaParameters.PAR2,
                         _OpcUaParameters.PAR3)
        rm0 = _readActs(_OpcUaParameters.RM1, _OpcUaParameters.RM2,
                        _OpcUaParameters.RM3)
        masked_ima0 = _interf.acq4d(n_frames, _ott)

        _ott.parab(pos_par + piston)
        # ott.refflat(pos_rm + piston)

        par1 = _readActs(_OpcUaParameters.PAR1, _OpcUaParameters.PAR2,
                         _OpcUaParameters.PAR3)
        rm1 = _readActs(_OpcUaParameters.RM1, _OpcUaParameters.RM2,
                        _OpcUaParameters.RM3)
        masked_ima1 = _interf.acq4d(n_frames, _ott)

        _ott.parab(pos_par - piston)
        # ott.refflat(pos_rm - piston)

        par2 = _readActs(_OpcUaParameters.PAR1, _OpcUaParameters.PAR2,
                         _OpcUaParameters.PAR3)
        rm2 = _readActs(_OpcUaParameters.RM1, _OpcUaParameters.RM2,
                        _OpcUaParameters.RM3)
        masked_ima2 = _interf.acq4d(n_frames, _ott)

        par = _np.array([par0, par1, par2])
        rm = _np.array([rm0, rm1, rm2])
        cubetto = _np.ma.dstack((masked_ima0, masked_ima1, masked_ima2))
        if cube is None:
            cube = cubetto
        else:
            cube = _np.ma.dstack((cube, cubetto))

        par_list.append(par)
        rm_list.append(rm)
        _pyfits.writeto(_os.path.join(dove, 'par.fits'),
                       _np.array(par_list), overwrite=True)
        _pyfits.writeto(_os.path.join(dove, 'rm.fits'),
                       _np.array(rm_list), overwrite=True)
        _pyfits.writeto(_os.path.join(dove, 'images.fits'),
                       cube.data, overwrite=True)
        _pyfits.append(_os.path.join(dove, 'images.fits'),
                      cube.mask.astype(int), overwrite=True)

    _ott.parab(pos_par)
    _ott.refflat(pos_rm)
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
    store_in_folder = _fold_name.CALIBRATION_ROOT_FOLDER
    save = _tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()
    par2rm = -2.05
    zern_vect = []
    parpos = []
    rmpos = []
    par0 = _ott.parab()
    rm0 = _ott.refflat()
    n2move = _np.array([3, 4])
    thedirection = _np.array([-1, 1])
    n_frames_alignment = 3
    tt_for_align = '20210111_152430'

    for k in n2move:
        for v in thedirection:

            for i in range(nstep):
                par1 = par0.copy()
                parmove = stepamp * i * v
                par1[k] += parmove
                print('Moving PAR[%d] by %d' % (k, parmove))
                _ott.parab(par1)
                rm1 = rm0.copy()
                rmmove = stepamp * i * v * par2rm
                rm1[k] += rmmove
                _ott.refflat(rm1)
                print('Moving RM[%d] by %d' % (k, rmmove))
                par_cmd, rm_cmd = _a.ott_alignment(n_frames_alignment, 1,
                                                  _np.array([0, 1]),
                                                  _np.array([3, 4]),
                                                  tt_for_align)
                par2 = _ott.parab()
                rm2 = _ott.refflat()
                masked_ima = _interf.acq4d(nframes, 0)
                name = _Timestamp.now() + '.fits'
                fits_file_name = _os.path.join(dove, name)
                _pyfits.writeto(fits_file_name, masked_ima.data)
                _pyfits.append(fits_file_name, masked_ima.mask.astype(int))

                coef, mat = _zernike.zernikeFit(masked_ima, _np.arange(10) + 1)
                zern_vect.append(coef)
                parpos.append(par2)
                rmpos.append(rm2)

                fits_file_name = _os.path.join(dove, 'zernike.fits')
                _pyfits.writeto(fits_file_name, _np.array(zern_vect),
                               overwrite=True)
                fits_file_name = _os.path.join(dove, 'PAR_positions.fits')
                _pyfits.writeto(fits_file_name, _np.array(parpos), overwrite=True)
                fits_file_name = _os.path.join(dove, 'RM_positions.fits')
                _pyfits.writeto(fits_file_name, _np.array(rmpos), overwrite=True)

    _ott.parab(par0)
    _ott.refflat(rm0)
    return tt


def _readActs(n1, n2, n3):
    act1 = _opc.get_position(n1)
    act2 = _opc.get_position(n2)
    act3 = _opc.get_position(n3)
    return _np.array([act1, act2, act3])


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
    hduList = _pyfits.open(deltapos_filepath)
    deltapos = hduList[0].data
    dx = deltapos[:, 0] * amp
    dy = deltapos[:, 1] * amp
    save = _tracking_number_folder.TtFolder(_fold_name.PISTON_TEST_ROOT_FOLDER)
    dove, tt = save._createFolderToStoreMeasurements()
    par0 = _ott.parab()
    rm0 = _ott.refflat()
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
            _ott.parab(par_new)
            _ott.refflat(rm_new)
            par_cmd, rm_cmd = _a.ott_alignment(n_frames_alignment, 1,
                                              _np.array([0, 1]),
                                              _np.array([3, 4]),
                                              tt_for_align)

        par = _ott.parab()
        rm = _ott.refflat()
        par_list.append(par)
        rm_list.append(rm)
        masked_ima0 = _interf.acq4d(n_frames_meas, _ott)
        par[2] += piston_value
        _ott.parab(par)
        masked_ima1 = _interf.acq4d(n_frames_meas, _ott)
        par[2] -= piston_value
        _ott.parab(par)
        diff = masked_ima1 - masked_ima0
        name = 'diff_%04d.fits' % i
        _interf.save_phasemap(dove, name, diff)
        coef, mat = _zernike.zernikeFit(diff, _np.arange(10) + 1)
        coef_list.append(coef)

        fits_file_name = _os.path.join(dove, 'Zernike.fits')
        _pyfits.writeto(fits_file_name, _np.array(coef_list), overwrite=True)
        fits_file_name = _os.path.join(dove, 'PAR_Positions.fits')
        _pyfits.writeto(fits_file_name, _np.array(par_list), overwrite=True)
        fits_file_name = _os.path.join(dove, 'RM_Positions.fits')
        _pyfits.writeto(fits_file_name, _np.array(rm_list), overwrite=True)
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
    image0 = _interf.acq4d(10, 0)
    delta_list = []
    sum_list = []
    coef_list = []
    for i in range(val_vec.size):
        _opc._setAct(act, val_vec[i])
        image = _interf.acq4d(10, 0)
        delta = image - image0
        delta_list.append(delta)
        coef, mat = _zernike.zernikeFit(delta, _np.arange(10) + 1)
        coef_list.append(coef)
        sum = _np.sqrt(coef[1] ** 2 + coef[2] ** 2)
        sum_list.append(sum)
    quad = _np.array(sum_list)
    _plt.plot(val_vec, quad, '-o')
    return quad, _np.array(coef_list), image0, delta_list


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
    par0 = _ott.parab()
    rm0 = _ott.refflat()
    delta_par = []
    delta_rm = []
    delta_object = []
    delta_object2 = []
    delta_object3 = []

    save = _tracking_number_folder.TtFolder(_fold_name.MAPPING_TEST_ROOT_FOLDER)
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

    file = open(_os.path.join(dove, 'MappingInfo.txt'), "a")
    file.write('Mapping object: ' + object_to_move + '\n')
    file.close()
    slide0 = _ott.slide()
    rslide0 = _ott.rslide()
    angle0 = _ott.angle()
    slide = -1
    rslide = -1
    angle = -1
    for i in range(n_iter):
        if shift[0] != 0:
            slide = _ott.slide(slide0 + ((i + 1) * shift[0]))
#             if slide==0:
#                 raise Exception('HOMING! PAR WIN')
        if shift[1] != 0:
            rslide = _ott.rslide(rslide0 + ((i + 1) * shift[1]))
#             if rslide==0:
#                 raise Exception('HOMING! RM WIN')
        if shift[2] != 0:
            angle = _ott.angle(angle0 + ((i + 1) * shift[2]))

        _time.sleep(5)
        par_cmd, rm_cmd = _a.ott_alignment(n_frames_alignment, 1,
                                          _np.array([0, 1]), _np.array([3, 4]),
                                          tt_for_align)
        image = _interf.acq4d(1)
        name = 'image_%04d.fits' % i
        _interf.save_phasemap(dove, name, image)
        par = _ott.parab()
        rm = _ott.refflat()
        delta_par.append(par - par0)
        delta_rm.append(rm - rm0)
        delta_object.append(slide - slide0)
        delta_object2.append(rslide - rslide0)
        delta_object3.append(angle - angle0)

        fits_file_name = _os.path.join(dove, 'delta_slide.fits')
        _pyfits.writeto(fits_file_name, _np.array(delta_object), overwrite=True)
        fits_file_name = _os.path.join(dove, 'delta_rslide.fits')
        _pyfits.writeto(fits_file_name, _np.array(delta_object2), overwrite=True)
        fits_file_name = _os.path.join(dove, 'delta_PAR_positions.fits')
        _pyfits.writeto(fits_file_name, _np.array(delta_par), overwrite=True)
        fits_file_name = _os.path.join(dove, 'delta_RM_positions.fits')
        _pyfits.writeto(fits_file_name, _np.array(delta_rm), overwrite=True)
        fits_file_name = _os.path.join(dove, 'delta_ANGLE_positions.fits')
        _pyfits.writeto(fits_file_name, _np.array(delta_object3), overwrite=True)

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
    _a.ott_alignment(n_images, 1, _np.array([0,1]), _np.array([3,4]), tt)
    par0 = _ott.parab()
    rm0 = _ott.refflat()
    # Set perturbation from perturbation_vec
    focus = perturbation_vec[0]
    tip = perturbation_vec[1]
    tilt = perturbation_vec[2]
    par1 = _np.copy(par0)
    par1[2] += focus
    par1[3] += tip
    par1[4] += tilt
    rm1 = _np.copy(rm0)
    rm1[3] += -tip*2.05
    rm1[4] += -tilt*2.05
    # Initialization
    zern_vec = _np.arange(1, 12, 1)
    coeff = []
    # Start image, perturbation and perturbed image
    image = _interf.acq4d(n_images)
    pippo = _zernike.zernikeFit(image, zern_vec)
    coeff.append(pippo[0])
    _ott.parab(par1)
    _ott.refflat(rm1)
    image = _interf.acq4d(n_images)
    pippo = _zernike.zernikeFit(image, zern_vec)
    coeff.append(pippo[0])
    # TipTilt pre-alignment
    if pre is True:
        _a.ott_alignment(n_images, 1, _np.array([0,1]), _np.array([3,4]), tt)
        image = _interf.acq4d(2)
        pippo = _zernike.zernikeFit(image, zern_vec)
        coeff.append(pippo[0])
    # First alignment
    _a.ott_alignment(n_images, 1, _np.array([0,1,2,3,4]),
                    _np.array([0,1,2,3,4]), tt)
    image = _interf.acq4d(2)
    pippo = _zernike.zernikeFit(image, zern_vec)
    coeff.append(pippo[0])
    # Second alignment
    _a.ott_alignment(n_images, 1, _np.array([0,1]), _np.array([3,4]), tt)
    image = _interf.acq4d(n_images)
    pippo = _zernike.zernikeFit(image, zern_vec)
    coeff.append(pippo[0])
    # Check image
    image = _interf.acq4d(n_images)
    pippo = _zernike.zernikeFit(image, zern_vec)
    coeff.append(pippo[0])

    coeff_matrix = _np.array(coeff)*1e9
    parend = _ott.parab()

    store_in_folder = _os.path.join(_fold_name.REPEATABILITY_ROOT_FOLDER,
                                   'Alignment')
    save = _tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()
    fits_file_name = _os.path.join(dove, 'zernike.fits')
    _pyfits.writeto(fits_file_name, coeff_matrix, overwrite=True)

    print(parend-par0)
    return coeff_matrix, tt

