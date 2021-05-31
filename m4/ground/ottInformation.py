'''
Authors
  - C. Selmi:  written in 2021
'''
import time
import os
import numpy as np
from m4.configuration.config import fold_name
from m4.configuration.ott_parameters import OpcUaParameters

class OttInformation():

    def __init__(self, ott, interf):
        """The constructor """
        self._ott = ott
        self._interf = interf
        self._angle = None
        self._rslide = None
        self._slide = None
        self._par = None
        self._rm = None
        self._temperature = None

    def acquireAndLogOttPositions(self):
        self._angle = self._ott.angleRotator.getPosition()
        self._rslide = self._ott.referenceMirrorSlider.getPosition()
        self._slide = self._ott.parabolaSlider.getPosition()
        self._par = self._ott.parabola.getPosition()
        self._rm = self._ott.referenceMirror.getPosition()
        self._temperature = self._ott.temperature.getTemperature()

        self._informationLog()

    def temperatureTimeHistory(self, n_measure, delay_in_seconds):
        temp_matrix = np.zeros((n_measure, OpcUaParameters.num_PT_sensor))
        for i in range(n_measure):
            temp_matrix[i, :] = self._ott.temperature.getTemperature()
            time.sleep(delay_in_seconds)
        return temp_matrix

    def interferogramTimeHistory(self, n_measure, delay_in_seconds):
        cube = None
        for i in range(n_measure):
            image = self._interf.acquire_phasemap()
            if cube is None:
                cube = image
            else:
                cube = np.ma.dstack((cube, image))
            time.sleep(delay_in_seconds)
        return cube

    def _informationLog(self):
        file_name = os.path.join(fold_name.LOG_ROOT_FOLDER,
                                 'OttInformationsLog.txt')
        file = open(file_name, 'a+')
        file.write('Angle = %f ' %self._angle)
        file.write('\n')
        file.write('RM slide = %f ' %self._rslide)
        file.write('\n')
        file.write('PAR slide = %f ' %self._slide)
        file.write('\n')
        file.write('RM position =')
        for i in range(self._rm.size):
            file.write(' %f ' %self._rm[i])
        file.write('\n')
        file.write('PAR position =')
        for i in range(self._par.size):
            file.write(' %f ' %self._par[i])
        file.write('\n')
        file.write('Temperature =')
        for i in range(self._temperature.size):
            file.write(' %f ' %self._temperature[i])
        file.write('\n')
        file.close()
