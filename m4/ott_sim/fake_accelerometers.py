'''
Authors
  - C. Selmi: written in 2020
'''
import logging
import os
import numpy as np
import h5py
from m4.configuration.config import fold_name
from m4.ground.timestamp import Timestamp
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
        for i in range(9):
            if i == 0:
                signal = np.column_stack((vector, vector))
            else:
                signal = np.column_stack((signal, vector))

        tt = Timestamp.now()
        name = tt + '.h5'
        final_destination = os.path.join(fold_name.ACC_ROOT_FOLDER, name)
        #simulatore non rebinnato
        hf = h5py.File(final_destination, 'w')
        hf.create_dataset('Accelerometers', data=signal[:,1:])
        hf.attrs['DT'] = OpcUaParameters.accelerometers_dt
        hf.attrs['ID'] = OpcUaParameters.accelerometers_plc_id
        hf.attrs['DIR'] = ['X', 'Z', 'Y', 'Z']
        hf.attrs['TIME'] = signal[:,0]
        hf.close()
        return name
