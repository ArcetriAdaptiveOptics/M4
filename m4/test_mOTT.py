'''
Author: C. Selmi
'''

import time
#import itertools as it
import numpy as np
from matplotlib import pyplot as plt
from m4.configuration.ott_parameters import OpcUaParameters
plt.figure(figsize=(5,5))

### TEST COMANDI PAR ###
def main(x0, y0, n_step, move):
    from m4.utils.opc_ua_controller import OpcUaController
    opc = OpcUaController()

#   x0 = 261
#   y0 = 68
    x, y = spiral(n_step, x0, y0)
    #plt.figure(figsize=(5,5))
    #plt.show()
    for i in range(x.size):
        par_position = readParPosition(opc)
        print(par_position)

        setParPosition(opc, 0, x[i], y[i])
        par_position = readParPosition(opc)
        print(par_position)

        if move==1:
            opc.move_object(OpcUaParameters.PAR_KIN)
        time.sleep(2)
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
    plt.plot(x,y,'-x', color='blue')


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

def spaz(opc, x0, y0, step, range):
    epsilon = 2 * range / step
    
    for i in range(step):
        y = (y0 - range) + (epsilon * step)

        setParPosition(opc, 0, x0-range, y)
        opc.move_object(OpcUaParameters.PAR_KIN)
        opc.wait_for_stop(OpcUaParameters.PAR_KIN)

        setParPosition(opc, 0, x0+range, y)
        opc.move_object(OpcUaParameters.PAR_KIN)
        opc.wait_for_stop(OpcUaParameters.PAR_KIN)



# def spiral(X, Y):
#     x = y = 0
#     dx = 0
#     dy = -1
#     x_list = []
#     y_list = []
#     for i in range(max(X, Y)**2):
#         if (-X/2 < x <= X/2) and (-Y/2 < y <= Y/2):
#             x_list.append(x)
#             y_list.append(y)
#         if x == y or (x < 0 and x == -y) or (x > 0 and x == 1-y):
#             dx, dy = -dy, dx
#             x, y = x+dx, y+dy
#     x = np.array(x_list)
#     y = np.array(y_list)
#     plt.figure(figsize=(5,5))
#     plt.plot(x, y, '.-')
#     return
