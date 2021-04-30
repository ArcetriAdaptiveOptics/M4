'''
Authors
  - C. Selmi: written in 2020
'''

import logging
import numpy as np
from m4.configuration.ott_parameters import OpcUaParameters
from m4.devices.base_temperature_sensors import BaseTemperatureSensors


class FakeTemperatureSensors(BaseTemperatureSensors):
    ''' Class for PT simulation
    '''

    def __init__(self):
        """The constructor """
        self._temp = np.zeros(OpcUaParameters.num_PT_sensor)
        self._logger = logging.getLogger('FakeParabola')

    def getTemperature(self):
        return self._temp
