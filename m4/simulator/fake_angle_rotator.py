"""
Authors
  - C. Selmi: written in 2020
"""

import logging


class FakeAngleRotator():
    """Class for ring angle rotation simulation (range: 0 to 360)

    HOW TO USE IT::

        from m4.ott_sim.fake_angle_rotator import FakeAngleRotator
        ang = FakeAngleRotator()
        angle = ang.getAngle()
        new_angle = ang.setAngle(absolute_position_in_deg)
    """

    def __init__(self):
        """The constructor"""
        self._angle = 0
        self._logger = logging.getLogger("FakeAngleRotator")

    def getPosition(self):
        """
        Returns
        -------
        angle: float
            angle position in degree
        """
        self._logger.debug("Position = %g" % self._angle)
        return self._angle

    def setPosition(self, absolute_position_in_deg):
        """
        Parameters
        ----------
        absolute_position_in_deg: float
            absolute position to set in degree
        """
        self._angle = absolute_position_in_deg
        return self.getPosition()
