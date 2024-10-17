'''
Authors
  - P. Ferraiuolo, R. Briguglio: written in 2024
'''
#import logging
import numpy as np
from m4.configuration.ott_parameters import OttParameters
from m4.devices.base_m4_exapode import BaseM4Exapode

'''
The kinematics is written in the form (normalization-less):

          C     ^Y
                |
       A     B  -->X
  Pist  Tx  Ty
A    1  -1  -1
B    1   1  -1
C    1   0   1
'''


class DpMotors(BaseM4Exapode):
    ''' Class for M4 control via opc ua??

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.m4_exapode import OpcUaM4Exapode
    '''
    DPinterface = 0 #dummy, to replicate the sintax 
    def __init__(self, DPinterface):
        """The constructor """
        self._dpinterface = DPinterface
        self._logger = logging.getLogger('DPMotors')
        self._kinematrix = np.array([[1,1,1],[-1,1,0],[-1,-1,1]])
        self._invkinematrix = np.linalg.inv(self._kinematrix)

    def getPosition(self):
        ''' to be implemented
        '''
        pp = _getMotorPosition()
        kk = _motor2kinematics(pp)
        pos = np.zeros(6)
        pos[OttParameters.M4_DOF] = kk[1:]
        print(pos)
        #return pos
        pass

    def setPosition(self, absolute_position_in_mm):
        ''' to be implemented
        '''
        v = np.array([0,absolute_position_in_mm[OttParameters.M4_DOF]])
        vm = np.dot(v,self._kinematrix)
        print(vm)
        _setMotorPosition(vm)
        pass


    def _getMotorPosition(self):
        ''' to be implemented, at low level
        '''
        #read_from_the_motors
        pass


    def _setMotorPosition(self, motorcmd):
        ''' to be implemented, at low level
        '''
        #command_the_motors
        pass

    def _motor2kinematics(self, pos):
        '''
        '''
        ptt = np.dot(pos,self._kinematrix)
        #return ptt
        pass

    def _kinematics2motor(self,pos):
        '''
        '''
        abc = np.dot(pos,self._invkinematrix)
        #return abc
        pass
    
    
