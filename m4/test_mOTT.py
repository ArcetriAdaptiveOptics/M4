'''
Author: C. Selmi
'''
import numpy as np
import itertools as it
from matplotlib import pyplot as plt
from m4.configuration.ott_parameters import OpcUaParameters


### TEST COMANDI PAR ###
def main():
    from m4.utils.opc_ua_controller import OpcUaController
    opc = OpcUaController()

    par_position = readParActsPosition(opc)
    print(par_position)

    setParActsPosition(opc, 0, 0, 0)
    par_position = readParActsPosition(opc)
    print(par_position)

    commandParActsPosition(opc)

    opc.stop_single_command(OpcUaParameters.PAR_KIN)


def readParActsPosition(opc):
    piston = opc.get_position(OpcUaParameters.PAR_PISTON)
    tip = opc.get_position(OpcUaParameters.PAR_TIP)
    tilt = opc.get_position(OpcUaParameters.PAR_TILT)
    return np.array([piston, tip, tilt])

def setParActsPosition(opc, piston, tip, tilt):
    opc.set_target_position(OpcUaParameters.PAR_PISTON, piston)
    opc.set_target_position(OpcUaParameters.PAR_TIP, tip)
    opc.set_target_position(OpcUaParameters.PAR_TILT, tilt)

def commandParActsPosition(opc):
    opc.move_object(OpcUaParameters.PAR_KIN)



def spirale(n):
    x = np.array([0])
    y = np.array([0])

    j = -1
    for k in range(0, n-1):
        for i in range(2):
            j += 1
            p0 = np.zeros(k+1)
            p1 = (p0+1)*(-1)**k
            p = np.array(list(it.accumulate(p1)))
            print(j)
            print(p0, p)
            if i == 0:
                x = np.append(x, p+x[j])
                y = np.append(y, p0+y[j])
            else:
                x = np.append(x, p0+x[j])
                y = np.append(y, p+y[j])

    plt.figure(figsize=(5,5))
