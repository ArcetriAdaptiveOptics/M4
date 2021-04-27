'''
Authors
  - C. Selmi: written in 2020
'''
import logging
from m4.devices.base_temperature_sensors import BaseTemperatureSensors

class OpcUaTemperatureSensors(BaseTemperatureSensors):
    ''' Class for Pt sensors control via opc ua
    '''

    def __init__(self, opcUa):
        """The constructor """
        self._opcUa = opcUa
        self._logger = logging.getLogger('OpcUaPT')

    def getTemperature(self):
        return self._opcUa.get_temperature_vector()
