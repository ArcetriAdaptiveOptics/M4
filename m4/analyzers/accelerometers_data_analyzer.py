'''
Authors
  - C. Selmi: written in 2021
'''

import numpy as np
from matplotlib import pyplot as plt
from m4.configuration.ott_parameters import OpcUaParameters
from m4.configuration import config_folder_names as fold_name
from m4.type.accelerometers_data import AccelerometersData


class AccelerometersDataAnalyzer():

    def __init__(self, tt):
        """The constructor """
        self.acc = AccelerometersData.loadInfoFromAccTtFolder(tt)
        self.tt = tt
        self._h5_file_path = self.acc._h5_file_path
        self.datah5 = None
        self._spe = None
        self._freq = None

    def readAndShow(self):
        self._spe, self._freq = self.acc.power_spectrum()
        self.datah5 = self.acc.datah5
        self.plot_power_spectrum()

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
