'''
Authors
  - C. Selmi: written in 2020
  '''
import logging
from m4.devices.base_parabola_slider import BaseParabolaSlider


class FakeParabolaSlider(BaseParabolaSlider):
    ''' Class for parabola slider simulation (range: -0.9 m +0.9 m)
    '''

    def __init__(self):
        """The constructor """
        self._pos = 0
        self._logger = logging.getLogger('FakeParabolaSlider')

    def getPosition(self):
        self._logger.debug('Position = %g' % self._pos)
        return self._pos

    def setPosition(self, absolute_position_in_mm):
        self._pos = absolute_position_in_mm
        return self.getPosition()

    def getPositionInM(self):
        return self._pos * 1e-3

    def setPositionInM(self, absolute_position_in_m):
        self._pos = absolute_position_in_m * 1e3
        return self.getPositionInM()
