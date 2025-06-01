"""
Authors
  - C. Selmi: written in 2020
"""

import logging


class OpcUaTemperatureSensors:
    """Class for Pt sensors control via opc ua

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.temperature_sensors import OpcUaTemperatureSensor
        sens = OpcUaTemperatureSensor(OpcUa)
    """

    def __init__(self, opcUa):
        """The constructor"""
        self._opcUa = opcUa
        self._logger = logging.getLogger("OpcUaPT")

    def getTemperature(self):
        """
        Returns
        -------
        temperature: numpy array
            vector containing temperature
        """
        return self._opcUa.get_temperature_vector()
