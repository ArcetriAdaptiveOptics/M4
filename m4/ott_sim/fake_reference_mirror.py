'''
Authors
  - C. Selmi: written in 2020
'''

import logging
import numpy as np
from m4.devices.base_reference_mirror import BaseReferenceMirror

class FakeReferenceMirror(BaseReferenceMirror):
    ''' Class for reference mirror simulation

    HOW TO USE IT::

        from m4.ott_sim.fake_reference_mirror import FakeReferenceMirror
        rm_slider = FakeReferenceMirror()
        pos = rm_slider.getPosition()
        new_pos = rm_slider.setPosition(absolute_position_in_mm)
    '''

    def __init__(self):
        """The constructor """
        self._pos = np.zeros(6)
        self._logger = logging.getLogger('FakeReferenceMirror')

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
            vector of six numbers containing dof values of rm

        Returns
        -------
        current_pos: numpy array [mm]
            absolute reference mirror position in millimeters
        '''
        self._pos = absolute_position_in_mm
        return self.getPosition()

    def getPositionInM(self):
        '''
        Returns
        -------
        current_pos: numpy array [m]
            reference mirror position in meters
        '''
        return self._pos*1e-3

    def setPositionInM(self, absolute_position_in_m):
        '''
        Parameters
        ----------
        absolute_position_in_mm: int [m]

        Returns
        -------
        current_pos: numpy array [m]
            absolute reference mirror position in meters
        '''
        self._pos = absolute_position_in_m * 1e3
        return self.getPositionInM()
