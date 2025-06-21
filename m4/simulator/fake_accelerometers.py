"""
Authors
  - C. Selmi: written in 2020
"""

import logging
import os
import numpy as np
import h5py
from m4.configuration import folders as fold_name
from opticalib.ground.osutils import newtn
from m4.configuration.ott_parameters import OpcUaParameters


class FakeAccelerometers():
    """Class for simulated accelerometers control

    HOW TO USE IT::

        from m4.ott_sim.fake_accelerometers import FakeAccelerometers
        acc = FakeAccelerometers()
        tt = acc. acquireData(recording_seconds)
    """

    def __init__(self):
        """The constructor"""
        self._logger = logging.getLogger("FakeAccelerometers")
        self._dt = OpcUaParameters.accelerometers_dt

    def acquireData(self, recording_seconds=5):
        """some function to simulate accelerometers data

        Parameters
        ----------
        recording_seconds: int [s]
            number of seconds for data recording

        Returns
        -------
        tt: string
            tracking number of mesurements
        """
        T = recording_seconds
        n = int(T / self._dt)
        t = np.linspace(0, T, n)
        freqSin = 2
        ampSin = 1
        vector = ampSin * np.sin(2 * np.pi * freqSin * t)
        for i in range(9):
            if i == 0:
                signal = np.column_stack((vector, vector))
            else:
                signal = np.column_stack((signal, vector))

        tt = newtc()
        name = tt + ".h5"
        final_destination = os.path.join(fold_name.ACCELEROMETERS_ROOT_FOLDER, name)
        # simulatore non rebinnato
        hf = h5py.File(final_destination, "w")
        hf.create_dataset("Accelerometers", data=signal[:, 1:])
        hf.attrs["DT"] = OpcUaParameters.accelerometers_dt
        hf.attrs["ID"] = OpcUaParameters.accelerometers_plc_id
        hf.attrs["DIR"] = ["X", "Z", "Y", "Z"]
        hf.attrs["TIME"] = signal[:, 0]
        hf.close()
        return tt
