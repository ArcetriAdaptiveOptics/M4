'''
Authors
  - C. Selmi: written in 2021
'''

import os
import numpy as np
import h5py
from m4.ground import rebinner
from m4.configuration.ott_parameters import OpcUaParameters
from m4.configuration import config_folder_names as fold_name


class AccelerometersData():
    '''
    '''

    def __init__(self):
        """The constructor """
        self._rebinnig_factor = OpcUaParameters.accelerometers_dt/OpcUaParameters.accelerometers_dt_plc
        self.tt = None
        self._h5_file_path = None
        self.datah5 = None
        self._spe = None
        self._freq = None

        self.dt = None
        self.id_vector = None
        self.directions = None
        self.time = None
        self.plc_voltscale = None
        self.plc_countscale = None
        self.sensitivity = None


# Functions for acquisition #
    def convertAndSaveData(self, start, final_destination):
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

    def counts_to_ms2(self, vec):
        id = OpcUaParameters.accelerometers_plc_id - 1
        sens = OpcUaParameters.accelerometers_sensitivity[id]
        plcfs = OpcUaParameters.accelerometers_plc_range[id]
        cal_list = []
        for i in range(vec.shape[1]):
            cal_vec = (vec[:, i]/OpcUaParameters.accelerometers_plc_totcounts)*plcfs*9.81/sens
            cal_list.append(cal_vec)
        cal_vec = np.array(cal_list)
        return cal_vec.T


    @staticmethod
    def loadInfoFromAccTtFolder(tt):
        """ Creates the object using information about path measurements

        Parameters
        ----------
                tt: string
                        measurement tracking number

        Returns
        -------
                theObject: object
                        analyzerIFF class object
        """
        theObject = AccelerometersData()
        theObject.tt = tt
        theObject._h5_file_path = os.path.join(
            fold_name.ACC_ROOT_FOLDER, tt + '.h5')
        hf = h5py.File(theObject._h5_file_path, 'r')
        try:
            theObject.dt = hf.attrs['DT']
            theObject.id_vector = hf.attrs['ID']
            theObject.directions = hf.attrs['DIR']
            theObject.time = hf.attrs['TIME']
            theObject.plc_voltscale = hf.attrs['PLC_VoltScale']
            theObject.plc_countscale = hf.attrs['PLC_CountScale']
            theObject.sensitivity = hf.attrs['Sensitivity']
        except Exception:
            theObject.dt = OpcUaParameters.accelerometers_dt_plc
            theObject.id_vector = OpcUaParameters.accelerometers_plc_id
            theObject.directions = ['X', 'Z', 'Y', 'Z']
        return theObject

# Functions for analysis #
    def read_data(self):
        '''
        Returns
        -------
            _datah5: numpy array
                measurements data
        '''
        if self._h5_file_path is None:
            raise OSError('Tracking number folder not specified ')
        hf = h5py.File(self._h5_file_path, 'r')
        self.datah5 = hf.get('Accelerometers')[()]
        return self.datah5

    def power_spectrum(self):
        '''
        Returns
        -------
            spe_list: list
                    list containing the spectrum of vectors composing
                    the matrix z
            freq_list: list
                    list containing the frequencies of vectors composing
                    the matrix z
        '''
        if self.datah5 is None:
            self.datah5 = self.read_data()

        if self.dt == OpcUaParameters.accelerometers_dt_plc:
            vec_cut = self.datah5[:, OpcUaParameters.accelerometers_plc_id]
            z = vec_cut.T
        else:
            z = self.datah5
        # dt = OpcUaParameters.accelerometers_dt
        # z = vec.T
    #         #spe = np.fft.fftshift(np.fft.rfft(vector, norm='ortho'))
    #         #freq = np.fft.fftshift(np.fft.rfftfreq(vector.size, d=self._dt))
        spe = np.fft.rfft(z, axis=1, norm='ortho')
        nn = np.sqrt(spe.shape[1])  # modRB
        self._spe = (np.abs(spe)) / nn
        self._freq = np.fft.rfftfreq(z.shape[1], d=self.dt)
        return self._spe, self._freq
