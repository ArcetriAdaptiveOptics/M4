'''
Authors
  - C. Selmi: written in 2020
'''
import logging
import numpy as np
from m4.configuration.ott_parameters import OpcUaParameters, OttParameters
from m4.devices.base_reference_mirror import BaseReferenceMirror

class OpcUaReferenceMirror(BaseReferenceMirror):
    ''' Class for reference mirror control via opc ua
    '''

    def __init__(self, opcUa):
        """The constructor """
        self._opcUa = opcUa
        self._logger = logging.getLogger('OpcUaReferenceMirror')

    def getPosition(self):
        current_pos = self._readRMPosition()
        self._logger.debug('Position = %s' % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_mm):
        n_opc = np.array([OpcUaParameters.RM_PISTON,
                          OpcUaParameters.RM_TIP,
                          OpcUaParameters.RM_TILT])
        for i in range(OttParameters.RM_DOF_PISTON.size):
            j = OttParameters.RM_DOF_PISTON[i]
            self._opcUa.set_target_position(n_opc[i], absolute_position_in_mm[j])
            # print(start_position[j])
        self._opcUa.move_object(OpcUaParameters.RM_KIN)
        self._opcUa.wait_for_stop(OpcUaParameters.RM_KIN)
        return self.getPosition()

    def _readRMPosition(self):
        piston = self._opcUa.get_position(OpcUaParameters.RM_PISTON)
        tip = self._opcUa.get_position(OpcUaParameters.RM_TIP)
        tilt = self._opcUa.get_position(OpcUaParameters.RM_TILT)
        return np.array([0, 0, piston, tip, tilt, 0])