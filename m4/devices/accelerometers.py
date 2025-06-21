"""
Authors
  - C. Selmi: written in 2020
"""

from opticalib.ground import logger as _l
import os as _os
import time as _time
import zmq as _zmq
from m4.configuration import folders as _fn
from opticalib.ground.osutils import newtn as _ntn
from m4.configuration.ott_parameters import OpcUaParameters as _opcua
from m4.type.accelerometers_data import AccelerometersData as _ad


class ZmqAccelerometers:
    """Class for accelerometers control via zmq

    HOW TO USE IT::

        from m4.devices.accelerometers import ZmqAccelerometers
        acc = ZmqAccelerometers()
        tt = acc. acquireData(recording_seconds)
    """

    def __init__(self):
        """The constructor"""
        self._logger = _l.set_up_logger("Accelerometers")
        self._acc = _ad()

    def acquireData(self, recording_seconds: int = 5):
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
        context = _zmq.Context()
        socket = context.socket(_zmq.REQ)
        socket.connect(_opcua.accelerometers_server)
        socket.send_string("START %d" % recording_seconds)
        _time.sleep(1)
        try:
            reply = socket.recv(1)
            print("Data from %s" % reply)
        except:
            raise OSError("No reply from socket")
        socket.disconnect(_opcua.accelerometers_server)

        lista = _os.listdir(_opcua.accelerometers_data_folder)
        lista.sort()
        h5_file_name = lista[len(lista) - 1]
        tt = _ntn()
        name = tt + ".h5"
        final_destination = _os.path.join(_fn.ACCELEROMETERS_ROOT_FOLDER, name)
        print("To %s" % final_destination)
        start = _os.path.join(_opcua.accelerometers_data_folder, h5_file_name)

        self._waitForEndAcquisition(start)
        self._acc.convertAndSaveData(start, final_destination)
        # self._tt = name
        return tt

    def _waitForEndAcquisition(self, data_file_path: str):
        t0 = _os.path.getmtime(data_file_path)
        _time.sleep(2)
        t1 = _os.path.getmtime(data_file_path)
        diff = t1 - t0
        while diff != 0:
            # print(diff)
            t0 = _os.path.getmtime(data_file_path)
            _time.sleep(2)
            t1 = _os.path.getmtime(data_file_path)
            diff = t1 - t0
