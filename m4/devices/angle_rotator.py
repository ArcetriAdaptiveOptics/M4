'''
Authors
  - C. Selmi: written in 2020
'''
import logging
from m4.configuration.ott_parameters import OpcUaParameters
from m4.devices.base_angle_rotator import BaseAngleRotator

class AngleRotator(BaseAngleRotator):
    ''' Class for ring angle rotation via opc ua '''

    def __init__(self, opcUa):
        """The constructor """
        self._opcUa = opcUa
        self._logger = logging.getLogger('OpcUaAngleRotator')

    def getPosition(self):
        current_pos = self._opcUa.get_position(OpcUaParameters.RA)
        self._logger.debug('Position = %g' % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_deg):
        self._checkAngle(absolute_position_in_deg)
        self._opcUa.set_target_position(OpcUaParameters.RA, absolute_position_in_deg)
        self._opcUa.move_object(OpcUaParameters.RA)
        self._opcUa.wait_for_stop(OpcUaParameters.RA)
        return self.getPosition()

    def _checkAngle(self, angle):
        if angle <= OpcUaParameters.min_angle or angle >= OpcUaParameters.max_angle:
            raise OSError(' The required angle is incorrect: %d' % angle)
