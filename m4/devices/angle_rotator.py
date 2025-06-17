"""
Author(s)
---------
- Chiara Selmi: written in 2020
- Pietro Ferraiuolo : refactored in 2025
"""

from opticalib.ground import logger as _l
from m4.configuration.ott_parameters import OpcUaParameters


class OpcUaAngleRotator:
    """Class for ring angle rotation control via opc ua

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.angle_rotator import OpcUaAngleRotator
        ang = OpcUaAngleRotator(OpcUa)
        pos = ang.getPosition()
        pos = ang.setPosition(absolute_position_in_deg)
    """

    def __init__(self, opcUa: object):
        """The constructor"""
        self._opcUa = opcUa
        self._logger = _l.set_up_logger("OpcUaAngleRotator")

    def getPosition(self):
        """Function to get the rotating ring angle

        Returns
        -------
            current_pos: int [deg]
                            rotating ring position
        """
        current_pos = self._opcUa.get_position(OpcUaParameters.RA)
        self._logger.debug("Position = %g" % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_deg: float):
        """Function to set the rotating ring angle (range: 0 to 360)

        Parameters
        ----------
            absolute_position_in_deg: int [deg]
                rotating ring position to set

        Returns
        -------
            current_pos: int [deg]
                        rotating ring position
        """
        self._checkAngle(absolute_position_in_deg)
        self._opcUa.set_target_position(OpcUaParameters.RA, absolute_position_in_deg)
        self._opcUa.move_object(OpcUaParameters.RA)
        self._opcUa.wait_for_stop(OpcUaParameters.RA)
        return self.getPosition()

    def _checkAngle(self, angle: float):
        """Function for input parameter control"""
        if angle <= OpcUaParameters.min_angle or angle >= OpcUaParameters.max_angle:
            raise OSError(" The required angle is incorrect: %d" % angle)
