"""
Authors
-------
Chiara Selmi: written in 2020
Pietro Ferraiuolo: updated in 2026 w/ new API
"""
import numpy as np
from opticalib.ground.logger import SystemLogger
from .opc_ua_controller import OpcUaController
from m4.configuration.ott_parameters import OpcUaParameters, OttParameters


class OpcUaParabola:
    """Class for parabola control via opc ua

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.parabola import OpcUaParabola
        par = OpcUaParabola(opcUa)
    """

    def __init__(self, opcUa: OpcUaController):
        """The constructor"""
        self._opcua = opcUa
        self._logger = SystemLogger(__class__)

    def getPosition(self):
        """
        Gets the parabola tripod's actuators position.

        Returns
        -------
        current_pos: int [-, -, mm, arcsec, arcsec, -]
            Tripod parabola position
        """
        ptt = self._opcua.get_act_positions(OpcUaParameters.PAR)
        current_pos = np.array(
            [0, 0, float(ptt[0]), float(ptt[1]), float(ptt[2]), 0]
        )
        self._logger.debug("Position = %s" % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_mm: list[float | int]) -> list[float]:
        """
        Function to set the absolute position of the parabola tripod's actuators.

        Parameters
        ----------
        absolute_position_in_mm: ArrayLike
            Vector of the six dof values of parabola.

        Returns
        -------
        print of the reached position
        """
        dofs = OttParameters.PARABOLA_DOF
        to_command = absolute_position_in_mm[dofs]
        self._opcua.set_act_positions(OpcUaParameters.PAR, to_command)
        self._opcua.wait_for_stop(OpcUaParameters.PAR)
        print(self.getPosition())
