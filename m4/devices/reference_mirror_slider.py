'''
Authors
  - C. Selmi: written in 2020
'''
from m4.configuration.ott_parameters import OpcUaParameters
import logging
from m4.devices.base_reference_mirror_slider import BaseReferenceMirrorSlider


class OpcUaReferenceMirrorSlider(BaseReferenceMirrorSlider):
    ''' Class for reference mirror slide control via opc ua
    '''

    def __init__(self, opcUa):
        """The constructor """
        self._opcUa = opcUa
        self._logger = logging.getLogger('OpcUaReferenceMirrorSlider')

    def getPosition(self):
        current_pos = self._opcUa.get_position(OpcUaParameters.CAR)
        self._logger.debug('Position = %s' % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_mm):
        self._checkRslide(absolute_position_in_mm)
        self._opcUa.set_target_position(OpcUaParameters.CAR,
                                        absolute_position_in_mm)
        self._opcUa.move_object(OpcUaParameters.CAR)
        self._opcUa.wait_for_stop(OpcUaParameters.CAR)
        return self.getPosition()

    def _checkRslide(self, r_slide):
        if r_slide <= OpcUaParameters.min_r_slide or r_slide >= OpcUaParameters.max_r_slide:
            raise OSError(' The required reference flat position is incorrect: %d' % r_slide)
