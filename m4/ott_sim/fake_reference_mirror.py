'''
Authors
  - C. Selmi: written in 2020
'''

import logging
import numpy as np
from m4.devices.base_reference_mirror import BaseReferenceMirror

class FakeReferenceMirror(BaseReferenceMirror):
    ''' Class for reference mirror simulation
    '''

    def __init__(self):
        """The constructor """
        self._pos = np.zeros(6)
        self._logger = logging.getLogger('FakeReferenceMirror')

    def getPosition(self):
        self._logger.debug('Position = %g' % self._pos)
        return self._pos

    def setPosition(self, absolute_position_in_mm):
        self._pos = absolute_position_in_mm
        return self.getPosition()
