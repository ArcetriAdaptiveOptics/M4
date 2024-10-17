'''
Authors
  - P. Ferraiuolo, R. Briguglio: written in 2024
'''
import logging
from m4.devices.base_m4_exapode import BaseM4Exapode

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

    def getPosition(self):
        ''' to be implemented
        '''
        pass

    def setPosition(self, absolute_position_in_mm):
        ''' to be implemented
        '''
        pass


    def getMotorPosition(self);
        ''' to be implemented, at low level
        '''
        pass


    def setMotorPosition(self, absolute_position_in_mm);
        ''' to be implemented, at low level
        '''
        pass

    def _dp_kinematics(self, pos):
        '''
        '''
        pass
    
    
