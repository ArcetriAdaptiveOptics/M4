'''
Authors
  - C. Selmi: written in 2020
'''
import logging
import os
import time
import zmq
from m4.configuration.ott_parameters import OpcUaParameters
from m4.devices.base_accelerometers import BaseAccelerometers

class ZmqAccelerometes(BaseAccelerometers):
    ''' Class for accelerometers control via zmq
    '''

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('Accelerometers')

    def acquireData(self, recording_seconds=5):
        '''
        Parameters
        ----------
        recording_seconds: int [s]
            recording seconds for data acquisition

        Returns
        -------
        measurement_data_folder: string
                        path of measurement data
        '''
        context = zmq.Context()
        socket = context.socket(zmq.REQ)
        socket.connect(OpcUaParameters.accelerometers_server)
        socket.send_string("START %d" %recording_seconds)
        time.sleep(1)
        try:
            reply = socket.recv(1)
            print('Data from %s' %reply)
        except:
            raise OSError('No reply from socket')
        socket.disconnect(OpcUaParameters.accelerometers_server)

        list = os.listdir(OpcUaParameters.accelerometers_data_folder)
        list.sort()
        h5_file_name = list[len(list)-1]
        measurement_data_folder = os.path.join(OpcUaParameters.accelerometers_data_folder,
                                               h5_file_name)
        return measurement_data_folder
