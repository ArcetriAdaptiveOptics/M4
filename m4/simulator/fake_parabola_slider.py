"""
Authors
  - C. Selmi: written in 2020
"""

import logging


class FakeParabolaSlider:
    """Class for parabola slider simulation (range: -0.9 m +0.9 m)

    HOW TO USE IT::

        from m4.ott_sim.fake_parabola_slider import FakeParabolaSlider
        par_slider = FakeParabolaSlider()
        pos = par_slider.getPosition()
        new_pos = par_slider.setPosition(absolute_position_in_mm)
    """

    def __init__(self):
        """The constructor"""
        self._pos = 0
        self._logger = logging.getLogger("FakeParabolaSlider")

    def getPosition(self):
        """
        Returns
        -------
        current_pos: int [mm]
            parabola slider position in millimeters
        """
        self._logger.debug("Position = %g" % self._pos)
        return self._pos

    def setPosition(self, absolute_position_in_mm):
        """
        Parameters
        ----------
        absolute_position_in_mm: int [mm]

        Returns
        -------
        current_pos: int [mm]
            absolute parabola slider position in millimeters
        """
        self._pos = absolute_position_in_mm
        return self.getPosition()

    def getPositionInM(self):
        """
        Returns
        -------
        current_pos: int [m]
            parabola slider position in meters
        """
        return self._pos * 1e-3

    def setPositionInM(self, absolute_position_in_m):
        """
        Parameters
        ----------
        absolute_position_in_mm: int [m]

        Returns
        -------
        current_pos: int [m]
            absolute parabola slider position in meters
        """
        self._pos = absolute_position_in_m * 1e3
        return self.getPositionInM()
