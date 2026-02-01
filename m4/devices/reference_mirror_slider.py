"""
Authors
  - C. Selmi: written in 2020
"""

from m4.configuration.ott_parameters import OpcUaParameters
from opticalib.ground.logger import SystemLogger
from .opc_ua_controller import OpcUaController

class OpcUaReferenceMirrorSlider:
    """
    Class for reference mirror slider control via opc ua

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.reference_mirror_slider import OpcUaReferenceMirrorSlider
        rm_slider = OpcUaReferenceMirrorSlider(opcUa)
    """

    def __init__(self, opcUa: OpcUaController):
        """The constructor"""
        self._opcUa = opcUa
        self._logger = SystemLogger(__class__)

    def getPosition(self) -> float:
        """
        Function to get the reference mirror slider position

        Returns
        -------
        current_pos: float
            reference mirror slider position, in millimeters
        """
        current_pos = self._opcUa.get_positions(OpcUaParameters.CAR)
        self._logger.debug("Position = %s" % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_mm: float):
        """
        Function to set the absolute position of the reference mirror slider

        Parameters
        ----------
        absolute_position_in_mm: float
            The absolute position of the reference mirror slider, in millimeters

        Returns
        -------
        current_pos: float
            Absolute reference mirror slider position
        """
        self._checkRslide(absolute_position_in_mm)
        self._opcUa.set_position(OpcUaParameters.CAR, absolute_position_in_mm)
        self._opcUa.wait_for_stop(OpcUaParameters.CAR)
        print(self.getPosition())

    def _checkRslide(self, r_slide: float):
        """Function for input parameter control"""
        if (
            r_slide <= OpcUaParameters.min_r_slide
            or r_slide >= OpcUaParameters.max_r_slide
        ):
            raise ValueError(
                " The required reference mirror slider position is incorrect: %d" % r_slide
            )

    # per OttImages.ottview
    def getPositionInM(self):
        """
        Returns
        -------
        current_pos: float
            reference mirror slider position in meters
        """
        pos = self.getPosition()
        return pos * 1e-3
