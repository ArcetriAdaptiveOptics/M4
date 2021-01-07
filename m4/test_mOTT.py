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
from m4.configuration import start
from m4.alignment import Alignment
from m4.ground.interface_4D import comm4d
from m4.ground import tracking_number_folder
from m4.ground import zernike
from m4.ground.timestamp import Timestamp

### TEST COMANDI PAR ###
def main(x0, y0, n_step, move):
    from m4.ground.opc_ua_controller import OpcUaController
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

def test_calib(commandAmpVector):
    ott = start.create_ott()
    a = Alignment(ott)

    nPushPull = 4
    n_frames = 5
    tt_list = []
    file_name = os.path.join(fold_name.CALIBRATION_ROOT_FOLDER, 'TtSleepCalib.txt')
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

def anlyzerTestCalib():
    tts1 = np.array(['20201214_091212', '20201214_092528', '20201214_093842',
                    '20201214_095152', '20201214_100508', '20201214_101821',
                    '20201214_103128', '20201214_104441', '20201214_105754',
                    '20201214_111110', '20201214_112435', '20201214_113749'])
    tts2 = np.array(['20201214_115451', '20201214_120323', '20201214_121200',
                     '20201214_122040', '20201214_122922', '20201214_123807',
                     '20201214_124640', '20201214_125504', '20201214_130327',
                     '20201214_131134', '20201214_131950', '20201214_132822'])
    intMat1 = None
    intMat2 = None
    for i in range(tts1.size):
        mat1 = pippo(tts1[i])
        mat2 = pippo(tts2[i])
        if intMat1 is None:
            intMat1 = mat1
        else:
            intMat1 = np.stack((intMat1, mat1))
        if intMat2 is None:
            intMat2 = mat2
        else:
            intMat2 = np.stack((intMat2, mat2))
    return intMat1, intMat2


def pippo(tt):
    from m4.utils.optical_alignment import opt_alignment
    al = opt_alignment(tt)
    intMat, rec, mask = al._loadAlignmentInfo()
    return intMat


def stability_test(n_images, delay):
    interf = comm4d()
    ott = start.create_ott()
    store_in_folder = fold_name.OPD_SERIES_ROOT_FOLDER
    save = tracking_number_folder.TtFolder(store_in_folder)
    dove, tt = save._createFolderToStoreMeasurements()

    vect_list = []
    t0 = time.time()
    for i in range(n_images):
        ti = time.time()
        dt = ti -t0
        masked_ima = interf.acq4d(ott, 1, 0)
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

def repeatability_test(n_meas, piston_value, n_frames):
    ott = start.create_ott()
    interf = comm4d()
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
        masked_ima0 = interf.acq4d(ott, n_frames, 0)

        ott.parab(pos_par + piston)
        #ott.refflat(pos_rm + piston)

        par1 = _readActs(OpcUaParameters.PAR1, OpcUaParameters.PAR2, OpcUaParameters.PAR3)
        rm1 = _readActs(OpcUaParameters.RM1, OpcUaParameters.RM2, OpcUaParameters.RM3)
        masked_ima1 = interf.acq4d(ott, n_frames, 0)

        ott.parab(pos_par - piston)
        #ott.refflat(pos_rm - piston)

        par2 = _readActs(OpcUaParameters.PAR1, OpcUaParameters.PAR2, OpcUaParameters.PAR3)
        rm2 = _readActs(OpcUaParameters.RM1, OpcUaParameters.RM2, OpcUaParameters.RM3)
        masked_ima2 = interf.acq4d(ott, n_frames, 0)

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

def readRepData(tt):
    file_name = os.path.join(fold_name.REPEATABILITY_ROOT_FOLDER, tt)
    hduList = pyfits.open(os.path.join(file_name, 'par.fits'))
    par = hduList[0].data
    hduList = pyfits.open(os.path.join(file_name, 'rm.fits'))
    rm = hduList[0].data
    #hduList = pyfits.open(os.path.join(file_name, 'images.fits'))
    #cube = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
    return par, rm

def analyzeOptRep(tt):
    par, rm, cube = readRepData(tt)
    z_list=[]
    for i in range(cube.shape[2]):
        masked_ima = cube[:,:,i]
        coef, mat = zernike.zernikeFit(masked_ima, np.arange(2, 7))
        z_list.append(coef)
    return np.array(z_list)

def analyserRepData(tt):
    par, rm = readRepData(tt)

    pos01_list_std = []
    pos02_list_std = []
    pos01_list_mean = []
    pos02_list_mean = []
    pos0_list = []
    for i in range(par.shape[2]):
        pos01 = par[:, 0, i] - par[:, 1, i]
        pos02 = par[:, 0, i] - par[:, 2, i]
        pos01_list_std.append(pos01.std())
        pos02_list_std.append(pos02.std())
        pos01_list_mean.append(pos01.mean())
        pos02_list_mean.append(pos02.mean())

        pos0 = par[:,0,i]
        pos0_list.append(pos0.std())

    pos01_std = np.array(pos01_list_std)
    pos02_std = np.array(pos02_list_std)
    pos01_mean = np.array(pos01_list_mean)
    pos02_mean = np.array(pos02_list_mean)
    pos0 = np.array(pos0_list)
    return pos01_std, pos02_std, pos01_mean, pos02_mean, pos0

def _readActs(n1, n2, n3):
    from m4.ground.opc_ua_controller import OpcUaController
    opc = OpcUaController()

    act1 = opc.get_position(n1)
    act2 = opc.get_position(n2)
    act3 = opc.get_position(n3)
    return np.array([act1, act2, act3])

def test_align():
    p0 = np.array([120, -114, -279,1319])
    p1 = np.array([194, 553, -438, -62])
    p = np.vstack((p0,p1))
    r0 = np.array([-9.2, -3])
    r1 = np.array([0.54, 5.8])
    r = np.vstack((r0,r1))

    ip = np.linalg.pinv(p)
    a = np.dot(ip, r)
    u,s,v = np.linalg.svd(a)
    return

def piston_test(piston_value, deltapos_filepath, amp, tt_for_align):
    # '/home/m4/pardeltapos.fits'
    hduList = pyfits.open(deltapos_filepath)
    deltapos = hduList[0].data
    dx = deltapos[:, 0] * amp
    dy = deltapos[:, 1] * amp
    ott = start.create_ott()
    interf = comm4d()
    a = Alignment(ott)
    save = tracking_number_folder.TtFolder(fold_name.PISTON_TEST_ROOT_FOLDER)
    dove, tt = save._createFolderToStoreMeasurements()
    par0 = ott.parab()
    rm0 = ott.refflat()
    n_frames = 5

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
            rm_new[3] += -2*dx[i]
            par_new[4] += dy[i]
            rm_new[4] += -2*dy[i]
            ott.parab(par_new)
            ott.refflat(rm_new)
            par_cmd, rm_cmd = a.ott_alignment(n_frames, 1, np.array([0,1]), np.array([3, 4]), tt_for_align)

        par = ott.parab()
        rm = ott.refflat()
        par_list.append(par)
        rm_list.append(rm)
        masked_ima0 = interf.acq4d(ott, n_frames)
        par[2] += piston_value
        ott.parab(par)
        masked_ima1 = interf.acq4d(ott, n_frames)
        par[2] -= piston_value
        ott.parab(par)
        diff = masked_ima1 - masked_ima0
        name = 'diff_%04d.fits' %i
        interf.save_phasemap(dove, name, diff)
        coef, mat = zernike.zernikeFit(diff, np.arange(3)+1)
        coef_list.append(coef)

        fits_file_name = os.path.join(dove, 'Zernike.fits')
        pyfits.writeto(fits_file_name, np.array(coef_list), overwrite=True)
        fits_file_name = os.path.join(dove, 'PAR_Positions.fits')
        pyfits.writeto(fits_file_name, np.array(par_list), overwrite=True)
        fits_file_name = os.path.join(dove, 'RM_Positions.fits')
        pyfits.writeto(fits_file_name, np.array(rm_list), overwrite=True)
        return tt
