'''
Authors
  - C. Selmi: written in 2021
'''

import os
import numpy as np
import h5py
from matplotlib import pyplot as plt
from m4.configuration.ott_parameters import OpcUaParameters
from m4.configuration.config import fold_name


class AccelerometersDataAnalyzer():

    def __init__(self, tt):
        """The constructor """
        self.tt = tt
        self._h5_file_path = os.path.join(
            AccelerometersDataAnalyzer._storageFolder(), tt + '.h5')
        hf = h5py.File(self._h5_file_path, 'r')
        try:
            self.dt = hf.attrs['DT']
            self.id_vector = hf.attrs['ID']
            self.directions = hf.attrs['DIR']
            self.time = hf.attrs['TIME']
            self.plc_voltscale = hf.attrs['PLC_VoltScale']
            self.plc_countscale = hf.attrs['PLC_CountScale']
            self.sensitivity = hf.attrs['Sensitivity']
        except Exception:
            self.dt = OpcUaParameters.accelerometers_dt_plc
            self.id_vector = OpcUaParameters.accelerometers_plc_id
            self.directions = ['X', 'Z', 'Y', 'Z']
        self._datah5 = None
        self._spe = None
        self._freq = None

    @staticmethod
    def _storageFolder():
        """ Creates the path where to save measurement data"""
        return fold_name.ACC_ROOT_FOLDER

    def readAndShow(self):
        spe, freq = self.power_spectrum()
        self.plot_power_spectrum()

    def _read_data(self):
        '''
        Returns
        -------
            _datah5: numpy array
                measurements data
        '''
        hf = h5py.File(self._h5_file_path, 'r')
        self._datah5 = hf.get('Accelerometers')[()]
        return self._datah5

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
        if self._datah5 is None:
            self._datah5 = self._read_data()

        if self.dt == OpcUaParameters.accelerometers_dt_plc:
            vec_cut = self._datah5[:, OpcUaParameters.accelerometers_plc_id]
            z = vec_cut.T
        else:
            z = self._datah5
        # dt = OpcUaParameters.accelerometers_dt
        # z = vec.T
    #         #spe = np.fft.fftshift(np.fft.rfft(vector, norm='ortho'))
    #         #freq = np.fft.fftshift(np.fft.rfftfreq(vector.size, d=self._dt))
        spe = np.fft.rfft(z, axis=1, norm='ortho')
        nn = np.sqrt(spe.shape[1])  # modRB
        self._spe = (np.abs(spe)) / nn
        self._freq = np.fft.rfftfreq(z.shape[1], d=self.dt)
        return self._spe, self._freq

    def plot_power_spectrum(self):
        spe1 = self._spe[:, 1:]
        freq1 = self._freq[1:]
        plt.figure()
        label_list = []

        for i in OpcUaParameters.accelerometers_plc_id:
            ss = 'Ch-' + str(i) + ' ' + OpcUaParameters.accelerometrs_directions[i - 1]
            label_list.append(ss)

        # label_list = OpcUaParameters.accelerometrs_directions[OpcUaParameters.accelerometers_plc_id]
        for i in range(spe1.shape[0]):
            plt.plot(freq1, np.abs(spe1[i,:]), '-', label=label_list[i])
        plt.xlabel('Freq[Hz]')
        plt.ylabel('Amplitude Spectrum |m/s2|')
        plt.xlim([0, 100])
        plt.title(self.tt)

        plt.ion()
        plt.show()
        plt.pause(0.01)
        plt.legend()
        plt.grid()
