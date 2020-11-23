'''
Author: C. Selmi
'''
import numpy as np
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
