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


    HOW TO USE IT::

        from m4.ott_sim.fake_temperature_sensor import FakeTemperatureSensor
        sens = FakeTemperatureSensor()
        temp = sens.getTemperature()
    '''

    def __init__(self):
        """The constructor """
        self._temp = np.zeros(OpcUaParameters.num_PT_sensor)
        self._logger = logging.getLogger('FakeParabola')

    def getTemperature(self):
        '''
        Returns
        -------
        temp: numpy array [C]
            vector cointaing temperature value of 20 sensors
        '''
        return self._temp
