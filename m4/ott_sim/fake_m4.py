'''
Authors
  - C. Selmi: written in 2020
'''

import numpy as np
import logging
from m4.devices.base_m4 import BaseM4

class FakeM4(BaseM4):
    ''' Class for M4 simulation
    '''

    def __init__(self):
        """The constructor """
        self._pos = np.zeros(6)
        self._logger = logging.getLogger('FakeM4')

    def getPosition(self):
        self._logger.debug('Position = %s' % self._pos)
        return self._pos

    def setPosition(self, absolute_position_in_mm):
        self._pos = absolute_position_in_mm
        return self.getPosition()

    def getPositionInM(self):
        return self._pos*1e-3

    def setPositionInM(self, absolute_position_in_m):
        self._pos = absolute_position_in_m * 1e3
        return self.getPositionInM()
