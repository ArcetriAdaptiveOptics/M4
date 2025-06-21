"""
Authors
  - C. Selmi:  written in 2021
"""

import time
import os
import numpy as np
from opticalib.ground.logger import txtLogger as _l
from m4.configuration import folders as _fn
from opticalib import typings as _ot
from m4.configuration.ott_parameters import OpcUaParameters


class OttInformation:
    """
    HOW TO USE IT::

        from m4.configuration import start
        ott, interf = start.create_ott('../M4/Data/SYSCONFData/Config.yaml')
        from m4.ground.ottInformation import OttInformation
        info = OttInformation(ott, interf)
    """

    def __init__(self, ott: _ot.GenericDevice, interf: _ot.InterferometerDevice):
        """The constructor"""
        self._ott = ott
        self._interf = interf
        self._angle = None
        self._rslide = None
        self._slide = None
        self._par = None
        self._rm = None
        self._temperature = None

    def acquireAndLogOttPositions(self):
        """Function for reading and logging of current test tower positions"""
        self._angle = self._ott.angleRotator.getPosition()
        self._rslide = self._ott.referenceMirrorSlider.getPosition()
        self._slide = self._ott.parabolaSlider.getPosition()
        self._par = self._ott.parabola.getPosition()
        self._rm = self._ott.referenceMirror.getPosition()
        self._temperature = self._ott.temperature.getTemperature()

        self._informationLog()

    def temperatureTimeHistory(self, n_measure: int, delay_in_seconds: int = 1) -> _ot.MatrixLike:
        """
        Parameters
        ----------
        n_measure: int
            number of measurements to acquire
        delay_in_seconds: int [s]
            time in seconds for delay between the measurements

        Returns
        -------
        temp_matrix: numpy array [n_measure, num_PT_sensor]
            meauserements matrix
        """
        temp_matrix = np.zeros((n_measure, OpcUaParameters.num_PT_sensor))
        for i in range(n_measure):
            temp_matrix[i, :] = self._ott.temperature.getTemperature()
            time.sleep(delay_in_seconds)
        return temp_matrix

    def interferogramTimeHistory(self, n_measure: int, delay_in_seconds: int) -> _ot.CubeData:
        """
        Parameters
        ----------
        n_measure: int
            number of measurements to acquire
        delay_in_seconds: int [s]
            time in seconds for delay between the measurements

        Returns
        -------
        cube: numpy masked array [pixel, pixel, n_measure]
            meauserements cube
        """
        cube_list = []
        for i in range(n_measure):
            image = self._interf.acquire_phasemap()
            cube_list.append(image)
            time.sleep(delay_in_seconds)
        cube = np.dstack(cube_list)
        return cube

    def _informationLog(self):
        """Function for logging"""
        logger = _l(_fn.LOGGING_ROOT_FOLDER, "OttInformationsLog.txt")
        logger.log("Angle = %f" % self._angle)
        logger.log("RM slide = %f" % self._rslide)
        logger.log("PAR slide = %f" % self._slide)
        logger.log("RM position = %s" % " ".join(f"{x:.2f}" for x in self._rm))
        logger.log("PAR position = %s" % " ".join(f"{x:.2f}" for x in self._par))
        logger.log("Temperature = %s" % " ".join(f"{x:.2f}" for x in self._temperature))
        # file = open(file_name, "a+")
        # file.write("Angle = %f " % self._angle)
        # file.write("\n")
        # file.write("RM slide = %f " % self._rslide)
        # file.write("\n")
        # file.write("PAR slide = %f " % self._slide)
        # file.write("\n")
        # file.write("RM position =")
        # for i in range(self._rm.size):
            # file.write(" %f " % self._rm[i])
        # file.write("\n")
        # file.write("PAR position =")
        # for i in range(self._par.size):
        #     file.write(" %f " % self._par[i])
        # file.write("\n")
        # file.write("Temperature =")
        # for i in range(self._temperature.size):
        #     file.write(" %f " % self._temperature[i])
        # file.write("\n")
        # file.close()
