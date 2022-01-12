'''
Authors
  - C. Selmi: written in 2020
'''

import logging
import numpy as np
from m4.devices.base_parabola import BaseParabola

class FakeParabola(BaseParabola):
    ''' Class for parabola simulation

    HOW TO USE IT::

        from m4.ott_sim.fake_parabola import FakeParabola
        par = FakeParabola()
        pos = par.getPosition()
        new_pos = par.setPosition(absolute_position_in_mm)
    '''

    def __init__(self):
        """The constructor """
        self._pos = np.zeros(6)
        self._logger = logging.getLogger('FakeParabola')

    def getPosition(self):
        '''
        Returns
        -------
        current_pos: numpy array [mm]
            parabola position in millimeters
        '''
        self._logger.debug('Position = %s' % self._pos)
        return self._pos

    def setPosition(self, absolute_position_in_mm):
        '''
        Parameters
        ----------
        absolute_position_in_mm: numpy array [mm]
            vector of six numbers containing dof values of parabola

        Returns
        -------
        current_pos: numpy array [mm]
            absolute parabola position in millimeters
        '''
        self._pos = absolute_position_in_mm
        return self.getPosition()

    def getPositionInM(self):
        '''
        Returns
        -------
        current_pos: numpy array [m]
            parabola position in meters
        '''
        return self._pos*1e-3

    def setPositionInM(self, absolute_position_in_m):
        '''
        Parameters
        ----------
        absolute_position_in_m: numpy array [m]
            vector of six numbers containing dof values of parabola

        Returns
        -------
        current_pos: numpy array [m]
            absolute parabola position in meters
        '''
        self._pos = absolute_position_in_m * 1e3
        return self.getPositionInM()
