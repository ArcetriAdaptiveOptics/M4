'''
Authors
  - C. Selmi: written in 2020
'''
import logging
from m4.devices.base_m4_exapode import BaseM4Exapode

class OpcUaM4Exapode(BaseM4Exapode):
    ''' Class for M4 control via opc ua??

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.m4_controller import OpcUaM4
    '''

    def __init__(self, opcUa):
        """The constructor """
        self._opcUa = opcUa
        self._logger = logging.getLogger('OpcUaM4')

    def getPosition(self):
        ''' to be implemented
        '''
        pass

    def setPosition(self, absolute_position_in_mm):
        ''' to be implemented
        '''
        pass
