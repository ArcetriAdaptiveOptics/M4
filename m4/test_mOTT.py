'''
Author: C. Selmi
'''

import time
import itertools as it
import numpy as np
from matplotlib import pyplot as plt
from m4.configuration.ott_parameters import OpcUaParameters


### TEST COMANDI PAR ###
def main():
    from m4.utils.opc_ua_controller import OpcUaController
    opc = OpcUaController()

    x0 = 0
    y0 = 1
    x, y = spiral(7, x0, y0)

    for i in range(x.size):
        par_position = readParPosition(opc)
        print(par_position)

        setParPosition(opc, 0, x[i], y[i])
        par_position = readParPosition(opc)
        print(par_position)

        opc.move_object(OpcUaParameters.PAR_KIN)
        time.sleep(2)



def readParPosition(opc):
    piston = opc.get_position(OpcUaParameters.PAR_PISTON)
    tip = opc.get_position(OpcUaParameters.PAR_TIP)
    tilt = opc.get_position(OpcUaParameters.PAR_TILT)
    return np.array([piston, tip, tilt])

def setParPosition(opc, piston, tip, tilt):
    opc.set_target_position(OpcUaParameters.PAR_PISTON, piston)
    opc.set_target_position(OpcUaParameters.PAR_TIP, tip)
    opc.set_target_position(OpcUaParameters.PAR_TILT, tilt)



def spiral(n, x0, y0):
    x = np.array([0])
    y = np.array([0])

    j = -1
    for k in range(0, n-1):
        for i in range(2):
            j = x.size - 1
            p0 = np.zeros(k+1)
            p1 = (p0+1)*(-1)**k
            p = np.array(list(it.accumulate(p1)))
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
    plt.figure(figsize=(5,5))
    plt.plot(x, y, '.-')
    return x, y

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
