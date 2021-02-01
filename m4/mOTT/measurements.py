'''
Authors
  - C. Selmi: written in 2020
'''

import time
import os
import numpy as np
from astropy.io import fits as pyfits
from matplotlib import pyplot as plt
from m4.configuration.ott_parameters import OpcUaParameters
from m4.ground.opc_ua_controller import OpcUaController
from m4.configuration.config import fold_name
from m4.configuration import start
from m4.alignment import Alignment
from m4.ground import tracking_number_folder
from m4.ground import zernike
from m4.ground.timestamp import Timestamp

ott, interf = start.create_ott()
a = Alignment(ott)
opc = OpcUaController()

def testCalib(commandAmpVector):
    nPushPull = 4
    n_frames = 5
    tt_list = []
    file_name = os.path.join(fold_name.CALIBRATION_ROOT_FOLDER,
                             'TtSleepCalib.txt')
    file = open(file_name, 'w+')
    for i in range(15):
        tt_tower = a.ott_calibration(n_frames, commandAmpVector, nPushPull, 0)
        tt_list.append(tt_tower)
        print(tt_tower)
        file.write('%s' %tt_tower)
    file.close()

    tt_list.append(commandAmpVector)
    ttAmpVector = np.array(tt_list)
    return ttAmpVector

def opticalMonitoring(n_images, delay):
    store_in_folder = fold_name.OPD_SERIES_ROOT_FOLDER
    save = tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()

    vect_list = []
    t0 = time.time()
    for i in range(n_images):
        ti = time.time()
        dt = ti -t0
        masked_ima = interf.acq4d(1, ott)
        name = Timestamp.now() + '.fits'
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, masked_ima.data)
        pyfits.append(fits_file_name, masked_ima.mask.astype(int))

        coef, mat = zernike.zernikeFit(masked_ima, np.arange(10)+1)
        vect = np.append(dt, coef)
        vect_list.append(vect)

        fits_file_name = os.path.join(dove, 'zernike.fits')
        pyfits.writeto(fits_file_name, np.array(vect_list), overwrite=True)

        time.sleep(delay)
    return tt

def actsRepeatability(n_meas, piston_value, n_frames):
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
        par0 = _readActs(OpcUaParameters.PAR1, OpcUaParameters.PAR2, OpcUaParameters.PAR3)
        rm0 = _readActs(OpcUaParameters.RM1, OpcUaParameters.RM2, OpcUaParameters.RM3)
        masked_ima0 = interf.acq4d(n_frames, ott)

        ott.parab(pos_par + piston)
        #ott.refflat(pos_rm + piston)

        par1 = _readActs(OpcUaParameters.PAR1, OpcUaParameters.PAR2, OpcUaParameters.PAR3)
        rm1 = _readActs(OpcUaParameters.RM1, OpcUaParameters.RM2, OpcUaParameters.RM3)
        masked_ima1 = interf.acq4d(n_frames, ott)

        ott.parab(pos_par - piston)
        #ott.refflat(pos_rm - piston)

        par2 = _readActs(OpcUaParameters.PAR1, OpcUaParameters.PAR2, OpcUaParameters.PAR3)
        rm2 = _readActs(OpcUaParameters.RM1, OpcUaParameters.RM2, OpcUaParameters.RM3)
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
        pyfits.writeto(os.path.join(dove, 'par.fits'), np.array(par_list), overwrite=True)
        pyfits.writeto(os.path.join(dove, 'rm.fits'), np.array(rm_list), overwrite=True)
        pyfits.writeto(os.path.join(dove, 'images.fits'), cube.data, overwrite=True)
        pyfits.append(os.path.join(dove, 'images.fits'), cube.mask.astype(int), overwrite=True)

    ott.parab(pos_par)
    ott.refflat(pos_rm)
    return tt

def scanAstigmComa(stepamp, nstep, nframes=10): #by RB 20210117. 
    #goal: to measure coma and astigmatism at different PAR position, spanning 500 arcsec
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
                parmove = stepamp*i*v
                par1[k] += parmove
                print('Moving PAR[%d] by %d' %(k, parmove))
                ott.parab(par1)
                rm1 = rm0.copy()
                rmmove = stepamp*i*v*par2rm
                rm1[k] += rmmove
                ott.refflat(rm1)
                print('Moving RM[%d] by %d' %(k, rmmove))
                par_cmd, rm_cmd = a.ott_alignment(n_frames_alignment, 1, np.array([0,1]), np.array([3, 4]), tt_for_align)
                par2 = ott.parab()
                rm2 = ott.refflat()
                masked_ima = interf.acq4d(nframes, 0)
                name = Timestamp.now() + '.fits'
                fits_file_name = os.path.join(dove, name)
                pyfits.writeto(fits_file_name, masked_ima.data)
                pyfits.append(fits_file_name, masked_ima.mask.astype(int))

                coef, mat = zernike.zernikeFit(masked_ima, np.arange(10)+1)
                zern_vect.append(coef)
                parpos.append(par2)
                rmpos.append(rm2)

                fits_file_name = os.path.join(dove, 'zernike.fits')
                pyfits.writeto(fits_file_name, np.array(zern_vect), overwrite=True)
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
            print('Iteration %d' %i)
            par_new = par0.copy()
            rm_new = rm0.copy()
            par_new[3] += dx[i]
            rm_new[3] += rmcoeff*dx[i]
            par_new[4] += dy[i]
            rm_new[4] += rmcoeff*dy[i]
            ott.parab(par_new)
            ott.refflat(rm_new)
            par_cmd, rm_cmd = a.ott_alignment(n_frames_alignment, 1, np.array([0,1]), np.array([3, 4]), tt_for_align)

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
        name = 'diff_%04d.fits' %i
        interf.save_phasemap(dove, name, diff)
        coef, mat = zernike.zernikeFit(diff, np.arange(10)+1)
        coef_list.append(coef)

        fits_file_name = os.path.join(dove, 'Zernike.fits')
        pyfits.writeto(fits_file_name, np.array(coef_list), overwrite=True)
        fits_file_name = os.path.join(dove, 'PAR_Positions.fits')
        pyfits.writeto(fits_file_name, np.array(par_list), overwrite=True)
        fits_file_name = os.path.join(dove, 'RM_Positions.fits')
        pyfits.writeto(fits_file_name, np.array(rm_list), overwrite=True)
    return tt

def parTiltTest(act, val_vec):
    image0 =interf.acq4d(10, 0)
    delta_list = []
    sum_list = []
    coef_list = []
    for i in range(val_vec.size):
        opc._setAct(act, val_vec[i])
        image = interf.acq4d(10, 0)
        delta = image - image0
        delta_list.append(delta)
        coef, mat = zernike.zernikeFit(delta, np.arange(10)+1)
        coef_list.append(coef)
        sum = np.sqrt(coef[1]**2+coef[2]**2)
        sum_list.append(sum)
    quad = np.array(sum_list)
    plt.plot(val_vec, quad, '-o')
    return quad, np.array(coef_list), image0, delta_list

def mappaturaPar(shift, n_iter, tt_for_align):
    '''
    shift = np.array([shift_par, shift_rm])
    '''
    n_object_to_move = OpcUaParameters.ST
    n_frames_alignment = 1
    par0 = ott.parab()
    rm0 = ott.refflat()
    delta_par = []
    delta_rm = []
    delta_object = []
    delta_object2 = []
    
    save = tracking_number_folder.TtFolder(fold_name.HOMING_TEST_ROOT_FOLDER)
    dove, tt = save._createFolderToStoreMeasurements()
    file = open(os.path.join(dove, 'HomingInfo.txt'), "a")
    file.write('Homing object number ' + str(n_object_to_move) +'\n')
    file.close()
    #n_iter = 20
    slide0 = ott.slide()
    rslide0 = ott.rslide()
    slide = -1
    rslide = -1
    for i in range(n_iter):
        if shift[0] != 0:
            slide = ott.slide(slide0+((i+1)*shift[0]))
            if slide==0:
                raise Exception('HOMING! PAR WIN')
        if shift[1] != 0:
            rslide = ott.rslide(rslide0+((i+1)*shift[1]))
            if rslide==0:
                raise Exception('HOMING! RM WIN')
        time.sleep(5)
        par_cmd, rm_cmd = a.ott_alignment(n_frames_alignment, 1,
                                          np.array([0,1]), np.array([3, 4]),
                                          tt_for_align)
        par = ott.parab()
        rm = ott.refflat()
        delta_par.append(par - par0)
        delta_rm.append(rm - rm0)
        delta_object.append(slide - slide0)
        delta_object2.append(rslide - rslide0)
        
        fits_file_name = os.path.join(dove, 'delta_slide.fits')
        pyfits.writeto(fits_file_name, np.array(delta_object), overwrite=True)
        fits_file_name = os.path.join(dove, 'delta_rslide.fits')
        pyfits.writeto(fits_file_name, np.array(delta_object2), overwrite=True)
        fits_file_name = os.path.join(dove, 'delta_PAR_positions.fits')
        pyfits.writeto(fits_file_name, np.array(delta_par), overwrite=True)
        fits_file_name = os.path.join(dove, 'delta_RM_positions.fits')
        pyfits.writeto(fits_file_name, np.array(delta_rm), overwrite=True)
    
    return tt     

         
     

                
                
    