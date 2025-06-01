"""
Authors
  - C. Selmi: written in 2020
"""

import logging
import numpy as np
from m4.configuration.ott_parameters import OpcUaParameters, OttParameters


class OpcUaParabola:
    """Class for parabola control via opc ua

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.parabola import OpcUaParabola
        par = OpcUaParabola(opcUa)
    """

    def __init__(self, opcUa):
        """The constructor"""
        self._opcUa = opcUa
        self._logger = logging.getLogger("OpcUaParabola")

    def getPosition(self):
        """Function to get the parabola position

        Returns
        -------
        current_pos: int [-, -, mm, arcsec, arcsec, -]
            parabola position
        """
        current_pos = self._readParPosition()
        self._logger.debug("Position = %s" % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_mm):
        """Function to set the absolute position of the parabola

        Parameters
        ----------
        absolute_position_in_mm: numpy array [-, -, mm, arcsec, arcsec, -]
            vector of six numbers containing dof values of parabola

        Returns
        -------
        current_pos: numpy array [mm]
            absolute parabola position
        """
        n_opc = np.array(
            [
                OpcUaParameters.PAR_PISTON,
                OpcUaParameters.PAR_TIP,
                OpcUaParameters.PAR_TILT,
            ]
        )
        for i in range(OttParameters.PARABOLA_DOF.size):
            j = OttParameters.PARABOLA_DOF[i]
            self._opcUa.set_target_position(n_opc[i], absolute_position_in_mm[j])
            # print(start_position[j])
        self._opcUa.move_object(OpcUaParameters.PAR_KIN)
        self._opcUa.wait_for_stop(OpcUaParameters.PAR_KIN)
        return self.getPosition()

    def _readParPosition(self):
        """
        Returns
        -------
            current_pos: numpy array [-, -, mm, arcsec, arcsec, -]
                        absolute parabola position
        """
        piston = self._opcUa.get_position(OpcUaParameters.PAR_PISTON)
        tip = self._opcUa.get_position(OpcUaParameters.PAR_TIP)
        tilt = self._opcUa.get_position(OpcUaParameters.PAR_TILT)
        return np.array([0, 0, float(piston), float(tip), float(tilt), 0])
