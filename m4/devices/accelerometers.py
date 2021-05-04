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
        self._rebinnig_factor = OpcUaParameters.accelerometers_dt/OpcUaParameters.accelerometers_dt_plc

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
        self._convertAndSaveData(start, final_destination)
        #self._tt = name
        return name

    def _convertAndSaveData(self, start, final_destination):
        hf = h5py.File(start, 'r')
        data = hf.get('Accelerometers')

        vec = data[:, OpcUaParameters.accelerometers_plc_id]
        #here starts modRB to implement data conversion all the time
        if self._rebinnig_factor != 1:
            rebinning_size = np.int(vec.shape[0]/self._rebinnig_factor)
            v_list=[]
            for i in range(vec.shape[1]):
                v_list.append(rebinner.rebin(vec[:, i], rebinning_size)) #1050

            out_vector = np.array(v_list)
            time = rebinner.rebin(data[:, 0], rebinning_size)

        else:
            out_vector = vec
            time = data[:, 0]

        nacc = out_vector.shape[0]
        print('Counts StDev:')
        for i in range(nacc):
            print(out_vector[i,:].std())

        vector_ms2 = self.counts_to_ms2(out_vector)

        hf = h5py.File(final_destination, 'w')
        hf.create_dataset('Accelerometers', data=vector_ms2)
        hf.attrs['DT'] = OpcUaParameters.accelerometers_dt
        hf.attrs['ID'] = OpcUaParameters.accelerometers_plc_id
        hf.attrs['DIR'] = OpcUaParameters.accelerometrs_directions
        hf.attrs['TIME'] = time
        hf.attrs['PLC_VoltScale'] = OpcUaParameters.accelerometers_plc_range
        hf.attrs['PLC_CountScale'] = OpcUaParameters.accelerometers_plc_totcounts
        hf.attrs['Sensitivity'] = OpcUaParameters.accelerometers_sensitivity
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

    def counts_to_ms2(self, vec):
        id = OpcUaParameters.accelerometers_plc_id - 1
        sens = OpcUaParameters.accelerometers_sensitivity[id]
        plcfs = OpcUaParameters.accelerometers_plc_range[id]
        cal_list = []
        for i in range(vec.shape[1]):
            cal_vec = (vec[:, i] /OpcUaParameters.accelerometers_plc_totcounts)*plcfs*9.81/sens
            cal_list.append(cal_vec)
        cal_vec = np.array(cal_list)
        return cal_vec.T
