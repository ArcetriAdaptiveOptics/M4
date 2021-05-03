'''
Authors
  - C. Selmi: written in 2020
'''
import logging
from m4.configuration.ott_parameters import OpcUaParameters
from m4.devices.base_angle_rotator import BaseAngleRotator

class OpcUaAngleRotator(BaseAngleRotator):
    ''' Class for ring angle rotation control via opc ua '''

    def __init__(self, opcUa):
        """The constructor """
        self._opcUa = opcUa
        self._logger = logging.getLogger('OpcUaAngleRotator')

    def getPosition(self):
        ''' Function to get the rotating ring angle

        Returns
        -------
            current_pos: int [deg]
                            rotating ring position
        '''
        current_pos = self._opcUa.get_position(OpcUaParameters.RA)
        self._logger.debug('Position = %g' % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_deg):
        ''' Function to set the rotating ring angle (range: 0 to 360)

        Parameters
        ----------
            absolute_position_in_deg: int [deg]
                rotating ring position to set

        Returns
        -------
            current_pos: int [deg]
                        rotating ring position
        '''
        self._checkAngle(absolute_position_in_deg)
        self._opcUa.set_target_position(OpcUaParameters.RA,
                                        absolute_position_in_deg)
        self._opcUa.move_object(OpcUaParameters.RA)
        self._opcUa.wait_for_stop(OpcUaParameters.RA)
        return self.getPosition()

    def _checkAngle(self, angle):
        ''' Function for input parameter control'''
        if angle <= OpcUaParameters.min_angle or angle >= OpcUaParameters.max_angle:
            raise OSError(' The required angle is incorrect: %d' % angle)
