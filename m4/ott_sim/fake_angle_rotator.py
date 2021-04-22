
'''
Authors
  - C. Selmi: written in 2020
'''
import logging
from m4.devices.base_angle_rotator import BaseAngleRotator


class FakeAngleRotator(BaseAngleRotator):
    ''' Class for ring angle rotation simulation
    '''

    def __init__(self):
        """The constructor """
        self._angle = 0
        self._logger = logging.getLogger('FakeAngleRotator')

    def getPosition(self):
        self._logger.debug('Position = %g' % self._angle)
        return self._angle

    def setPosition(self, absolute_position_in_deg):
        self._angle = absolute_position_in_deg
        return self.getPosition()
