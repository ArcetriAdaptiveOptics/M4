'''
Authors
  - C. Selmi: written in 2020
'''
import logging
import os
import time
import numpy as np
import h5py
import zmq
from m4.ground import rebinner
from m4.configuration.config import fold_name
from m4.ground.timestamp import Timestamp
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
        name: string
            tracking number of measurements
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
        tt = Timestamp.now()
        name = tt + '.h5'
        final_destination = os.path.join(fold_name.ACC_ROOT_FOLDER, name)
        print( 'To %s' %final_destination)
        start = os.path.join(OpcUaParameters.accelerometers_data_folder,
                             h5_file_name)

        self._waitForEndAcquisition(start)

        rebinnig_factor = OpcUaParameters.accelerometers_dt/OpcUaParameters.accelerometers_dt_plc
        if rebinnig_factor == 1:
            os.system('cp %s %s' %(start, final_destination)) #popen
        else:
            self._rebinAndSaveData(start, final_destination, rebinnig_factor)
        self._tt = name
        return name

    def _rebinAndSaveData(self, start, final_destination, rebinnig_factor):
        hf = h5py.File(start, 'r')
        data = hf.get('Accelerometers')

        vec = data[:, OpcUaParameters.accelerometers_plc_id]
        rebinning_size = np.int(vec.shape[0]/rebinnig_factor)
        v_list=[]
        for i in range(vec.shape[1]):
            v_list.append(rebinner.rebin(vec[:, i], rebinning_size)) #1050
        rebinned_vector = np.array(v_list)
        time_vector = rebinner.rebin(data[:, 0], rebinning_size)

        hf = h5py.File(final_destination, 'w')
        hf.create_dataset('Accelerometers', data=rebinned_vector)
        hf.attrs['DT'] = OpcUaParameters.accelerometers_dt
        hf.attrs['ID'] = OpcUaParameters.accelerometers_plc_id
        hf.attrs['DIR'] = ['X', 'Z', 'Y', 'Z']
        hf.attrs['TIME'] = time_vector
        hf.close()

    def _waitForEndAcquisition(self, data_file_path):
        t0 = os.path.getmtime(data_file_path)
        time.sleep(2)
        t1 = os.path.getmtime(data_file_path)
        diff = t1-t0
        while diff != 0:
            #print(diff)
            t0 = os.path.getmtime(data_file_path)
            time.sleep(2)
            t1 = os.path.getmtime(data_file_path)
            diff = t1-t0
