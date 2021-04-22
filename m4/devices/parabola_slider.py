'''
Authors
  - C. Selmi: written in 2020
'''
import logging
from m4.configuration.ott_parameters import OpcUaParameters
from m4.devices.base_parabola_slider import BaseParabolaSlider


class OpcUaParabolaSlider(BaseParabolaSlider):
    ''' Class for parabola slide control via opc ua
    '''

    def __init__(self, opcUa):
        """The constructor """
        self._opcUa = opcUa
        self._logger = logging.getLogger('OpcUaParabolaSlider')

    def getPosition(self):
        current_pos = self._opcUa.get_position(OpcUaParameters.ST)
        self._logger.debug('Position = %g' % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_mm):
        self._checkSlide(absolute_position_in_mm)
        self._opcUa.set_target_position(
            OpcUaParameters.ST, absolute_position_in_mm)
        self._opcUa.move_object(OpcUaParameters.ST)
        self._opcUa.wait_for_stop(OpcUaParameters.ST)
        return self.getPosition()

    def _checkSlide(self, slide):
        if slide <= OpcUaParameters.min_slide or slide >= OpcUaParameters.max_slide:
            raise ValueError(
                'The required parabola position is incorrect: %g' % slide)
