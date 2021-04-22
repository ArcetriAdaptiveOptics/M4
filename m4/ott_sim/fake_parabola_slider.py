import logging
from m4.devices.base_parabola_slider import BaseParabolaSlider


class FakeParabolaSlider(BaseParabolaSlider):

    def __init__(self):
        self._pos = 0
        self._logger = logging.getLogger('FakeParabolaSlider')

    def getPosition(self):
        self._logger.debug('Position = %g' % self._pos)
        return self._pos

    def setPosition(self, absolute_position_in_mm):
        self._pos = absolute_position_in_mm
        return self.getPosition()
