'''
Authors
  - C. Selmi: written in 2020
'''

import logging
import numpy as np
from m4.devices.base_parabola import BaseParabola

class FakeParabola(BaseParabola):
    ''' Class for parabola simulation
    '''

    def __init__(self):
        """The constructor """
        self._pos = np.zeros(6)
        self._logger = logging.getLogger('FakeParabola')

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
