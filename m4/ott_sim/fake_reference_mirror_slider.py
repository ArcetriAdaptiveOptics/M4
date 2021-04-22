'''
Authors
  - C. Selmi: written in 2020
'''
import logging
from m4.devices.base_reference_mirror_slider import BaseReferenceMirrorSlider

class FakeReferenceMirrorSlider(BaseReferenceMirrorSlider):
    ''' Class for reference mirror slider simulation
    '''

    def __init__(self):
        """The constructor """
        self._pos = 0
        self._logger = logging.getLogger('FakeReferenceMirrorSlider')

    def getPosition(self):
        self._logger.debug('Position = %g' % self._pos)
        return self._pos

    def setPosition(self, absolute_position_in_mm):
        self._pos = absolute_position_in_mm
        return self.getPosition()
