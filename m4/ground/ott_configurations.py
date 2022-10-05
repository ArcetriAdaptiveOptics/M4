'''
Authors
  - C. Selmi:  written in 2022
'''

class OttConfigurations():
    '''
    HOW TO USE IT::

        from m4.configuration import start
        ott, interf = start.create_ott('../M4/Data/SYSCONFData/Config.yaml')
        from m4.ground.ott_configurations import OttConfigurations
        oc = OttConfigurations(ott)
    '''

    def __init__(self, ott):
        """The constructor """
        self._ott = ott
        self._angle = None
        self._rslide = None
        self._pslide = None


    def move_to_segment_view(self, number_of_segment, RM_in):
        '''move the ott configuration to a specific section of the DM

        Parameters
        ----------
            number_of_segment: integer
                number of the target section
            RM: boolean
                Reference mirror in (True) or not (False)
        '''
        self._ott.parabolaSlider.setPosition(844)
        if RM_in==True:
            self._ott.referenceMirrorSlider.setPosition(844)
        if RM_in==False:
            self._ott.referenceMirrorSlider.setPosition(0)

        self._ott.angleRotator.setPosition(30+60*(number_of_segment-1))
        return

    def move_to_central_view(self, RM_in):
        self._ott.parabolaSlider.setPosition(0)
        if RM_in==True:
            self._ott.referenceMirrorSlider.setPosition(0)
        if RM_in==False:
            self._ott.referenceMirrorSlider.setPosition(999)

        self._ott.angleRotator.setPosition(0)
        return

    def get_configuration(self):
        '''
        Returns
        -------
        segment_view = boolean
            in the ott is in segment view configuration it is True,
            else False
        RM_in = boolean
            if reference mirror is inside the image it is True,
            else False
        '''
        self._angle = self._ott.angleRotator.getPosition()
        self._rslide = self._ott.referenceMirrorSlider.getPosition()
        self._pslide = self._ott.parabolaSlider.getPosition()

        if self._pslide>800:
            print('segment view')
            segment_view = True
            if self._rslide>800:
                print('reference mirror in')
                rm_in = True
            elif self._rslide==0:
                print('reference mirror out')
                rm_in = False
        if self._pslide==0:
            print('central view')
            segment_view = False
            if self._rslide==0:
                print('reference mirror in')
                rm_in = True
            elif self._rslide>=999:
                print('reference mirror out')
                rm_in = False
        return segment_view, rm_in
