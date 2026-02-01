"""
Authors
  - C. Selmi: written in 2020
"""
from opticalib.ground.logger import SystemLogger
from .opc_ua_controller import OpcUaController
from m4.configuration.ott_parameters import OpcUaParameters


class OpcUaParabolaSlider:
    """Class for parabola slide control via opc ua

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.parabola_slider import OpcUaParabolaSlider
        par_slider = OpcUaParabolaSlider(opcUa)
    """

    def __init__(self, opcUa: OpcUaController):
        """The constructor"""
        self._opcUa = opcUa
        self._logger = SystemLogger(__class__)

    def getPosition(self) -> float :
        """Function to get the parabola slider position

        Returns
        -------
        current_pos: int [mm]
            parabola slider position
        """
        current_pos = self._opcUa.get_positions(OpcUaParameters.ST)
        self._logger.debug("Position = %g" % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_mm: int):
        """
        Function to set the absolute position of the parabola slider

        Parameters
        ----------
        absolute_position_in_mm: int
            The position to reach, in millimeters

        Returns
        -------
        Prints the reached position
        """
        self._checkSlide(absolute_position_in_mm)
        self._opcUa.set_position(OpcUaParameters.ST, absolute_position_in_mm)
        self._opcUa.wait_for_stop(OpcUaParameters.ST)
        print(self.getPosition())

    def _checkSlide(self, slide: int):
        """Function for input parameter control"""
        if slide <= OpcUaParameters.min_slide or slide >= OpcUaParameters.max_slide:
            raise ValueError(
                "The required parabola slider position is incorrect: %g" % slide
            )

    # per OttImages.ottview
    def getPositionInM(self) -> float:
        """
        Returns
        -------
        current_pos: int [m]
            parabola slider position in meters
        """
        pos: float = self.getPosition()
        return pos * 1e-3
