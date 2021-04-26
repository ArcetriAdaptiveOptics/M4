'''
Authors
  - C. Selmi: written in 2020
'''

import numpy as np
import logging
from m4.devices.base_parabola import BaseParabola

class FakeParabolaSlider(BaseParabola):
    ''' Class for parabola simulation
    '''

    def __init__(self):
        """The constructor """
        self._pos = np.zeros(6)
        self._logger = logging.getLogger('FakeParabola')

    def getPosition(self):
        self._logger.debug('Position = %g' % self._pos)
        return self._pos

    def setPosition(self, absolute_position_in_mm):
        self._pos = absolute_position_in_mm
        return self.getPosition()
