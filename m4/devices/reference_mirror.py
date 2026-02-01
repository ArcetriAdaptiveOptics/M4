"""
Authors
  - C. Selmi: written in 2021
"""

from opticalib.ground.logger import SystemLogger
import numpy as np
from .opc_ua_controller import OpcUaController
from m4.configuration.ott_parameters import OpcUaParameters, OttParameters


class OpcUaReferenceMirror:
    """Class for reference mirror control via opc ua


    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.reference_mirror import OpcUaReferenceMirror
        rm = OpcUaReferenceMirror(opcUa)
    """

    def __init__(self, opcUa: OpcUaController):
        """The constructor"""
        self._opcUa = opcUa
        self._logger = SystemLogger(__class__)

    def getPosition(self):
        """
        Gets the reference mirror's actuators position.

        Returns
        -------
        current_pos: int [-, -, mm, arcsec, arcsec, -]
            Reference mirror position
        """
        ptt = self._opcUa.get_act_positions(OpcUaParameters.RM)
        current_pos = np.array(
            [0, 0, float(ptt[0]), float(ptt[1]), float(ptt[2]), 0]
        )
        self._logger.debug("Position = %s" % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_mm: list[float | int]) -> list[float]:
        """
        Function to set the absolute position of the reference mirror's actuators.

        Parameters
        ----------
        absolute_position_in_mm: ArrayLike
            Vector of the six dof values of reference mirror.

        Returns
        -------
        print of the reached position
        """
        dofs = OttParameters.RM_DOF_PISTON
        to_command = absolute_position_in_mm[dofs]
        self._opcUa.set_act_positions(OpcUaParameters.RM, to_command)
        self._opcUa.wait_for_stop(OpcUaParameters.RM)
        print(self.getPosition())

