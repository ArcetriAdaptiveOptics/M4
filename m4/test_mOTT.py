'''
Author: C. Selmi
'''

import time
import os
#import itertools as it
import numpy as np
from astropy.io import fits as pyfits
from IPython.display import clear_output
from matplotlib import pyplot as plt
from m4.configuration.ott_parameters import OpcUaParameters
from m4.configuration.config import fold_name
#plt.figure(figsize=(5,5))
from m4.configuration import start
from m4.alignment import Alignment
from m4.utils.interface_4D import comm4d
from m4.ground import tracking_number_folder
from m4.utils import image_extender as ie
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from m4.ground.timestamp import Timestamp

### TEST COMANDI PAR ###
def main(x0, y0, n_step, move):
    from m4.utils.opc_ua_controller import OpcUaController
    opc = OpcUaController()

#   coordinate
    x, y = spiral(n_step, x0, y0)
#     no_move = 'non muovo'
#     x, y = spaz(opc, x0, y0, step, intervallo, no_move)
    plt.figure(figsize=(5, 5))
    #plt.show()
    for i in range(x.size):
        par_position = readParPosition(opc)
        print(par_position)

        setParPosition(opc, 0, x[i], y[i])
        par_position = readParPosition(opc)
        print(par_position)

        if move==1:
            opc.move_object(OpcUaParameters.PAR_KIN)
            opc.wait_for_stop(OpcUaParameters.PAR_KIN)
        time.sleep(1)
        plotthespiral(x[0:i], y[0:i])



def readParPosition(opc):
    piston = opc.get_position(OpcUaParameters.PAR_PISTON)
    tip = opc.get_position(OpcUaParameters.PAR_TIP)
    tilt = opc.get_position(OpcUaParameters.PAR_TILT)
    return np.array([piston, tip, tilt])

def setParPosition(opc, piston, tip, tilt):
    opc.set_target_position(OpcUaParameters.PAR_PISTON, piston)
    opc.set_target_position(OpcUaParameters.PAR_TIP, tip)
    opc.set_target_position(OpcUaParameters.PAR_TILT, tilt)

def plotthespiral(x,y):
    clear_output(wait=True)
    plt.plot(x,y,'-x', color='blue')
    plt.show()
    plt.pause(0.1)


def spiral(n, x0, y0):
    x = np.array([0])
    y = np.array([0])

    j = -1
    for k in range(0, n-1):
        for i in range(2):
            j = x.size - 1
            p0 = np.zeros(k+1)
            p1 = (p0+1)*(-1)**k
            p = np.cumsum(p1)
#             print(j)
#             print(p0, p)
            if i == 0:
                x = np.append(x, p+x[j])
                y = np.append(y, p0+y[j])
            else:
                x = np.append(x, p0+x[j])
                y = np.append(y, p+y[j])

    x = x + x0
    y = y + y0
    #plt.figure(figsize=(5,5))
    #plt.plot(x, y, '.-')
    return x, y

def spaz(opc, x0, y0, step, intervallo, move=None):
    epsilon = 2 * intervallo / step

    x_list = []
    y_list = []
    for i in range(0, step):
        y = (y0 - intervallo) + (epsilon * i)

        if move is None:
            setParPosition(opc, 0, x0-intervallo, y)
            opc.move_object(OpcUaParameters.PAR_KIN)
            opc.wait_for_stop(OpcUaParameters.PAR_KIN)
        x_list.append(x0-intervallo)
        y_list.append(y)

        if move is None:
            setParPosition(opc, 0, x0+intervallo, y)
            opc.move_object(OpcUaParameters.PAR_KIN)
            opc.wait_for_stop(OpcUaParameters.PAR_KIN)
        x_list.append(x0+intervallo)
        y_list.append(y)
    return np.array(x_list), np.array(y_list)

### TEST ###

def test_calib(commandAmpVector, nPushPull, move):
    ott = start.create_ott()
    a = Alignment(ott)

    if move == 1:
        tt_tower = a.ott_calibration(commandAmpVector, nPushPull, 0)
#         tt_tower = a._cal.measureCalibrationMatrix(ott, 0, commandAmpVector,
#                                                       nPushPull)
        return tt_tower
    else:
        mat = a._cal._createCommandMatrix(0, commandAmpVector, nPushPull)
        plt.clf()
        plt.imshow(mat, origin='lower')
        plt.colorbar()
        return mat

def pippo(tt):
    from m4.utils.optical_alignment import opt_alignment
    al = opt_alignment(tt)
    intMat, rec, mask = al._loadAlignmentInfo()
    return intMat, rec, mask


def stability_test(n_images, delay):
    interf = comm4d()
    ott = start.create_ott()
    zOnM4 = ZernikeOnM4()
    store_in_folder = fold_name.OPD_SERIES_ROOT_FOLDER
    save = tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()

    vect_list = []
    t0 = time.time()
    for i in range(n_images):
    	ti = time.time()
    	dt = ti -t0
        masked_ima = interf.acq4d(ott, 1, 0)
        name = Timestamp.now()
        fits_file_name = os.path.join(dove, name)
        pyfits.writeto(fits_file_name, masked_ima.data)
        pyfits.append(fits_file_name, masked_ima.mask.astype(int))

        new_ima = ie.imageExtender(masked_ima)
        coef, mat = zOnM4.zernikeFit(new_ima, np.arange(2, 12))
        vect = np.append(dt, coef)
        vect_list.append(vect)
        
        fits_file_name = os.path.join(dove, 'zernike.fits')
        pyfits.writeto(fits_file_name, np.array(vect_list), overwrite=True)

        time.sleep(delay)

    return tt

def repeatability_test(n_meas, piston_value):
    ott = start.create_ott()
    interf = comm4d()
    store_in_folder = fold_name.REPEATABILITY_ROOT_FOLDER
    save = tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()

    piston = np.array([0, 0, piston_value, 0, 0, 0])

    for i in range(n_meas):
        par0 = _readActs(OpcUaParameters.PAR1, OpcUaParameters.PAR2, OpcUaParameters.PAR3)
        rm0 = _readActs(OpcUaParameters.RM1, OpcUaParameters.RM2, OpcUaParameters.RM3)
        masked_ima0 = interf.acq4d(ott, 1, 0)

        ott.parab(piston)
        ott.refflat(piston)

        par1 = _readActs(OpcUaParameters.PAR1, OpcUaParameters.PAR2, OpcUaParameters.PAR3)
        rm1 = _readActs(OpcUaParameters.RM1, OpcUaParameters.RM2, OpcUaParameters.RM3)
        masked_ima1 = interf.acq4d(ott, 1, 0)

        ott.parab(-piston)
        ott.refflat(-piston)

        par2 = _readActs(OpcUaParameters.PAR1, OpcUaParameters.PAR2, OpcUaParameters.PAR3)
        rm2 = _readActs(OpcUaParameters.RM1, OpcUaParameters.RM2, OpcUaParameters.RM3)
        masked_ima2 = interf.acq4d(ott, 1, 0)

        par = np.array([par0, par1, par2])
        rm = np.array([rm0, rm1, rm2])
        cube = np.ma.dstack((masked_ima0, masked_ima1, masked_ima2))
        pyfits.writeto(os.path.join(dove, 'par_%d.fits' %i), par)
        pyfits.writeto(os.path.join(dove, 'rm_%d.fits' %i), rm)
        pyfits.writeto(os.path.join(dove, 'images_%d.fits' %i), cube.data)
        pyfits.append(os.path.join(dove, 'images_%d.fits' %i), cube.mask.astype(int))
    return tt


def _readActs(n1, n2, n3):
    from m4.utils.opc_ua_controller import OpcUaController
    opc = OpcUaController()
    
    act1 = opc.get_position(n1)
    act2 = opc.get_position(n2)
    act3 = opc.get_position(n3)
    return np.array([act1, act2, act3])
