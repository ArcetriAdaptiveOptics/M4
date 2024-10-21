"""
Author(s)
---------
    - Pietro Ferraiuolo: written in 2024
    - Runa Briguglio: written in 2024
"""
import zmq
import numpy as np
from m4.ground import logger_set_up as logger
from m4.configuration.ott_parameters import OttParameters
from m4.devices.base_m4_exapode import BaseM4Exapode

# The kinematics is written in the form (normalization-less):

#           C     ^Y
#                 |
#        A     B  -->X
#   Pist  Tx  Ty
# A    1  -1  -1
# B    1   1  -1
# C    1   0   1

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
        self._kinematrix = np.array([[1,1,1],[-1,1,0],[-1,-1,1]])
        self._invkinematrix = np.linalg.inv(self._kinematrix)
        self._context = None
        self._socket = None
        self.remote_ip = 0
        self.remote_port = 0
        logger.set_up_logger('path/dpMotors.log', 10)

    def getPosition(self):
        ''' to be implemented
        '''
        pp = self._getMotorPosition()
        kk = self._motor2kinematics(pp)
        pos = np.zeros(6)
        pos[OttParameters.M4_DOF] = kk[1:]
        print(pos)
        #return pos
        pass

    def setPosition(self, absolute_position_in_mm):
        ''' to be implemented
        '''
        v = np.array([0,absolute_position_in_mm[OttParameters.M4_DOF]])
        vm = np.dot(v, self._kinematrix)
        print(vm)
        self._setMotorPosition(vm)
        pass


    def _getMotorPosition(self):
        ''' to be implemented, at low level
        '''
        self._connectBusBox()
        # Si manda il comando
        # si riceve risposta 
        # si fanno le conversioni del caso
        self._disconnectBusBox()
        pass


    def _setMotorPosition(self, motorcmd):
        ''' to be implemented, at low level
        '''
        self._connectBusBox()
        # Si manda il comando
        # si riceve risposta 
        # si fanno le conversioni del caso
        self._disconnectBusBox()
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
    
    def _connectBusBox(self):
        try:
            self._context = zmq.Context()
            self._socket = self._context.socket(zmq.PAIR)
            self._socketconnect(f"tcp://{self.remote_ip}:{self.remote_port}")
            return True
        except Exception as e:
            logger.log(f"ConnectZMQ: {e}")
            return False

    def _disconnectBusBox(self):
        try:
            self._socket.disconnect(f"tcp://{self.remote_ip}:{self.remote_port}")
            self._socket.close()
            self._context.term()
            return True
        except Exception as e:
            logger.log(f"DisconnectZMQ: {e}")
            return False
        pass
