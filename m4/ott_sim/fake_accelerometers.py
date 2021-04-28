'''
Authors
  - C. Selmi: written in 2020
'''
import logging
import numpy as np
from m4.configuration.ott_parameters import OpcUaParameters
from m4.devices.base_accelerometers import BaseAccelerometers

class FakeAccelerometers(BaseAccelerometers):
    ''' Class for simulated accelerometers control
    '''

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('FakeAccelerometers')
        self._dt = OpcUaParameters.accelerometers_dt

    def acquireData(self, recording_seconds=5):
        ''' some function to simulate accelerometers data '''
        T = recording_seconds
        n = int(T/self._dt)
        t = np.linspace(0, T, n)
        freqSin = 2
        ampSin = 1
        vector = ampSin * np.sin(2*np.pi*freqSin*t)
        for i in range(7):
            if i == 0:
                signal = np.column_stack((vector, vector))
            else:
                signal = np.column_stack((signal, vector))
        return signal
