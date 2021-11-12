'''
Authors
  - C. Selmi: written in 2021
'''

import numpy as np
from matplotlib import pyplot as plt
from m4.configuration.ott_parameters import OpcUaParameters
from m4.type.accelerometers_data import AccelerometersData


class AccelerometersDataAnalyzer():
    ''' Class used to analyze accelerometer data
    '''

    def __init__(self, tt):
        """The constructor """
        self._acc = AccelerometersData.loadInfoFromAccTtFolder(tt)
        self.tt = tt
        self._h5_file_path = self._acc._h5_file_path
        self.datah5 = None
        self._spe = None
        self._freq = None

    def readAndShow(self):
        '''
        Function for calculation and displaying power spectrum from accelerometers data
        '''
        self._spe, self._freq = self._acc.power_spectrum()
        self.datah5 = self._acc.datah5
        self.plot_power_spectrum()

    def getData(self):
        if self.datah5 is None:
            self.datah5 = self._acc.datah5
        return self.datah5

    def getSpecAndFreq(self):
        if self._spe is None and self._freq is None:
            self._spe, self._freq = self._acc.power_spectrum()
        return self._spe, self._freq

    def plot_power_spectrum(self):
        '''
        Function for displaying power spectrum
        '''
        spe, freq = self.getSpecAndFreq()
        spe1 = spe[:, 1:]
        freq1 = freq[1:]
        plt.figure()
        label_list = []

        for i in OpcUaParameters.accelerometers_plc_id:
            ss = 'Ch-' + str(i) + ' ' + OpcUaParameters.accelerometrs_directions[i - 1]
            label_list.append(ss)

# label_list = OpcUaParameters.accelerometrs_directions[OpcUaParameters.accelerometers_plc_id]
        for i in range(spe1.shape[0]):
            plt.plot(freq1, np.abs(spe1[i, :]), '-', label=label_list[i])
        plt.xlabel('Freq[Hz]')
        plt.ylabel('Amplitude Spectrum |m/s2|')
        plt.xlim([0, 100])
        plt.title(self.tt)

        plt.ion()
        plt.show()
        plt.pause(0.01)
        plt.legend()
        plt.grid()
