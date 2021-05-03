'''
Authors
  - C. Selmi: written in 2020
'''
import logging
from m4.devices.base_m4 import BaseM4

class OpcUaM4(BaseM4):
    ''' Class for M4 control via opc ua??
    '''

    def __init__(self, opcUa):
        """The constructor """
        self._opcUa = opcUa
        self._logger = logging.getLogger('OpcUaM4')

    def getPosition(self):
        pass

    def setPosition(self, absolute_position_in_mm):
        pass
