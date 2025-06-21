"""
Authors
  - C. Selmi: written in 2020
"""

import logging


class FakeReferenceMirrorSlider:
    """Class for reference mirror slider simulation (range: -0.05 m to 0.4 m)

    HOW TO USE IT::

        from m4.ott_sim.fake_reference_mirror_slider import FakereferenceMirrorSlider
        rm_slider = FakeReferenceMirrorSlider()
        pos = rm_slider.getPosition()
        new_pos = rm_slider.setPosition(absolute_position_in_mm)
    """

    def __init__(self):
        """The constructor"""
        self._pos = 0
        self._logger = logging.getLogger("FakeReferenceMirrorSlider")

    def getPosition(self):
        """
        Returns
        -------
        current_pos: int [mm]
            reference mirror slider position in millimeters
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
            absolute reference mirror slider position in millimeters
        """
        self._pos = absolute_position_in_mm
        return self.getPosition()

    def getPositionInM(self):
        """
        Returns
        -------
        current_pos: int [m]
            reference mirror slider position in meters
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
            absolute reference mirror slider position in meters
        """
        self._pos = absolute_position_in_m * 1e3
        return self.getPositionInM()
