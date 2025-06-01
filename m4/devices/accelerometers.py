"""
Authors
  - C. Selmi: written in 2020
"""

import logging
import os
import time
import zmq
from m4.configuration import config_folder_names as fold_name
from m4.configuration.ott_parameters import OpcUaParameters
from m4.type.accelerometers_data import AccelerometersData


class ZmqAccelerometers:
    """Class for accelerometers control via zmq

    HOW TO USE IT::

        from m4.devices.accelerometers import ZmqAccelerometers
        acc = ZmqAccelerometers()
        tt = acc. acquireData(recording_seconds)
    """

    def __init__(self):
        """The constructor"""
        self._logger = logging.getLogger("Accelerometers")
        self._acc = AccelerometersData()

    def acquireData(self, recording_seconds=5):
        """
        Parameters
        ----------
        recording_seconds: int [s]
            recording seconds for data acquisition

        Returns
        -------
        name: string
            tracking number of measurements
        """
        context = zmq.Context()
        socket = context.socket(zmq.REQ)
        socket.connect(OpcUaParameters.accelerometers_server)
        socket.send_string("START %d" % recording_seconds)
        time.sleep(1)
        try:
            reply = socket.recv(1)
            print("Data from %s" % reply)
        except:
            raise OSError("No reply from socket")
        socket.disconnect(OpcUaParameters.accelerometers_server)

        lista = os.listdir(OpcUaParameters.accelerometers_data_folder)
        lista.sort()
        h5_file_name = lista[len(lista) - 1]
        tt = Timestamp.now()
        name = tt + ".h5"
        final_destination = os.path.join(fold_name.ACC_ROOT_FOLDER, name)
        print("To %s" % final_destination)
        start = os.path.join(OpcUaParameters.accelerometers_data_folder, h5_file_name)

        self._waitForEndAcquisition(start)
        self._acc.convertAndSaveData(start, final_destination)
        # self._tt = name
        return tt

    def _waitForEndAcquisition(self, data_file_path):
        t0 = os.path.getmtime(data_file_path)
        time.sleep(2)
        t1 = os.path.getmtime(data_file_path)
        diff = t1 - t0
        while diff != 0:
            # print(diff)
            t0 = os.path.getmtime(data_file_path)
            time.sleep(2)
            t1 = os.path.getmtime(data_file_path)
            diff = t1 - t0
