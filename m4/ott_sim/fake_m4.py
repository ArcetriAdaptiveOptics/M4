'''
Authors
  - C. Selmi: written in 2020
'''

import numpy as np
import logging
from m4.devices.base_m4 import BaseM4

class FakeM4(BaseM4):
    ''' Class for M4 simulation

    HOW TO USE IT::

        from m4.ott_sim.fake_m4 import FakeM4
        dm = FakeM4()
        dm_pos = dm.getPosition()
        dm_pos = dm.setPosition()
    '''

    def __init__(self):
        """The constructor """
        self._pos = np.zeros(6)
        self._logger = logging.getLogger('FakeM4')

    def getPosition(self):
        '''
        Returns
        ------
        pos: nump array [mm]
            vector containing dof position of deformable mirror in millimiters
        '''
        self._logger.debug('Position = %s' % self._pos)
        return self._pos

    def setPosition(self, absolute_position_in_mm):
        '''
        Parameters
        ----------
        absolute_position_in_mm: numpy arrauy [mm]
            absolute position of m4 dof to set in millimeters
        '''
        self._pos = absolute_position_in_mm
        return self.getPosition()

    def getPositionInM(self):
        '''
        Returns
        ------
        pos: nump array [m]
            vector containing dof position of deformable mirror in meters
        '''
        return self._pos*1e-3

    def setPositionInM(self, absolute_position_in_m):
        '''
        Parameters
        ----------
        absolute_position_in_m: numpy arrauy [m]
            absolute position of m4 dof to set in meters
        '''
        self._pos = absolute_position_in_m * 1e3
        return self.getPositionInM()
