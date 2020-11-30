'''
Author: C. Selmi
'''

import time
#import itertools as it
import numpy as np
from IPython.display import clear_output
from matplotlib import pyplot as plt
from m4.configuration.ott_parameters import OpcUaParameters
#plt.figure(figsize=(5,5))
from m4.configuration import start
from m4.alignment import Alignment

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

def test_calib(commandAmpVector, nPushPull=None, move=None):
    ott = start.create_ott()
    a = Alignment(ott)
    if nPushPull is None:
        nPushPull = 3
    if move is None:
        tt_tower = a.ott_calibration(commandAmpVector, nPushPull, 3)
        return tt_tower
    else:
        mat = a._cal._createCommandMatrix(0, commandAmpVector, nPushPull)
        plt.clf()
        plt.imshow(mat, origin='lower')
        plt.colorbar()
        return mat
