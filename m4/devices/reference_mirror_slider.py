'''
Authors
  - C. Selmi: written in 2020
'''
import os
from m4.configuration.ott_parameters import OpcUaParameters
import logging
from m4.devices.base_reference_mirror_slider import BaseReferenceMirrorSlider
import playsound
from m4.configuration.ott_parameters import Sound

class OpcUaReferenceMirrorSlider(BaseReferenceMirrorSlider):
    ''' Class for reference mirror slider control via opc ua

    HOW TO USE IT::

        from m4.devices.opc_ua_controller import OpcUaController
        opcUa = OpcUaController()
        from m4.devices.reference_mirror_slider import OpcUaReferenceMirrorSlider
        rm_slider = OpcUaReferenceMirrorSlider(opcUa)
    '''

    def __init__(self, opcUa):
        """The constructor """
        self._opcUa = opcUa
        self._logger = logging.getLogger('OpcUaReferenceMirrorSlider')

    def getPosition(self):
        ''' Function to get the reference mirror slider position

        Returns
        -------
            current_pos: int [mm]
                            reference mirror slider position
        '''
        current_pos = self._opcUa.get_position(OpcUaParameters.CAR)
        self._logger.debug('Position = %s' % current_pos)
        return current_pos

    def setPosition(self, absolute_position_in_mm):
        '''Function to set the absolute position of the reference mirror slider

        Parameters
        ----------
        absolute_position_in_mm: int [mm]

        Returns
        -------
            current_pos: int [mm]
                        absolute reference mirror slider position
        '''
        self._checkRslide(absolute_position_in_mm)
        self._opcUa.set_target_position(OpcUaParameters.CAR,
                                        absolute_position_in_mm)
        self._opcUa.move_object(OpcUaParameters.CAR)
        self._opcUa.wait_for_stop(OpcUaParameters.CAR)
        if Sound.PLAY is True:
            playsound.playsound(os.path.join(Sound.AUDIO_FILE_PATH, 'ref-in.mp3'))
        return self.getPosition()

    def _checkRslide(self, r_slide):
        ''' Function for input parameter control'''
        if r_slide <= OpcUaParameters.min_r_slide or r_slide >= OpcUaParameters.max_r_slide:
            raise OSError(' The required reference flat position is incorrect: %d' % r_slide)
    
    #per OttImages.ottview
    def getPositionInM(self):
        '''
        Returns
        -------
        current_pos: int [m]
            reference mirror slider position in meters
        '''
        pos = self.getPosition()
        return pos*1e-3
