'''
Authors
  - C. Selmi: written in 2020
'''

class OTT():

    def __init__(self, parabola_slider, reference_mirror_slider, angle_rotator,
                 parabola, reference_mirror, m4, temperature_sensor):
        """The constructor """
        self._parabola_slider = parabola_slider
        self._reference_mirror_slider = reference_mirror_slider
        self._angle_rotator = angle_rotator
        self._parabola = parabola
        self._reference_mirror = reference_mirror
        self._m4 = m4
        self._pt = temperature_sensor


# Elements position
    @property
    def slide(self):
        ''' Function to set the parabola translation (range: -0.9 m +0.9 m)

        Other Parameters
        ----------------
        par_trans: int, optional [mm]
                If par_trans is not set return slider position

        Returns
        -------
            par_trans: int [mm]
                    parabola translation
        '''
        return self._parabola_slider.getPosition()
    @slide.setter
    def slide(self, par_slider_position):
        self._parabola_slider.setPosition(par_slider_position)

    def rslide(self, ref_flat_slider_position=None):
        '''  Function to set the reference flat mirror (range: -0.05 m to 0.4 m)

        Other Parameters
        ----------
        ref_flat: int, optional [mm]
                If ref_flat is not set return slider position

        Returns
        -------
        ref_flat: int [mm]
                reference flat mirror slider position
        '''
        if ref_flat_slider_position is None:
            return self._reference_mirror_slider.getPosition()
        else:
            return self._reference_mirror_slider.setPosition(ref_flat_slider_position)

    def angle(self, rot_ring_angle=None):
        ''' Function to set the rotating ring angle (range: 0 to 360)

        Other Parameters
        ---------------
            rot_ring_angle: int, optional [deg]
                If rot_ring_angle is not set returns rotating ring position

        Returns
        -------
            rot_ring_angle: int [deg]
                            rotating ring position
        '''
        if rot_ring_angle is None:
            return self._angle_rotator.getPosition()
        else:
            return self._angle_rotator.setPosition(rot_ring_angle)


# Elements alignment
    def parab(self, par_position=None):
        '''Function to set the absolute position of the parable

        Other Parameters
        ----------
        start_position: numpy array, optional [mm]
                        vector of six position,
                        If start_position is not set returns parabola position

        Returns
        -------
            start_position: numpy array [mm]
                        absolute parabola position

        '''
        if par_position is None:
            return self._parabola.getPosition()
        else:
            return self._parabola.setPosition(par_position)

    def refflat(self, rm_position=None):
        '''Function to set the absolute position of the reference flat

        Other Parameters
        ----------
        start_position: numpy array, optional [mm]
                        vector of six position,
                        If start_position is not set returns reference flat position

        Returns
        -------
            start_position: numpy array
                        absolute position of the reference flat
        '''
        if rm_position is None:
            return self._reference_mirror.getPosition()
        else:
            return self._reference_mirror.setPosition(rm_position)

    def m4(self, m4_position=None):
        '''Function to set deformable mirror's position

        Other Parameters
        ----------------
        start_position: numpy array, optional
                        vector of six position,
                        If start_position is not set it's equal to zero vector

        Returns
        -------
            start_position: numpy array
                        start position of the deformable mirror
        '''
        if m4_position is None:
            return self._m4.getPosition()
        else:
            return self._m4.setPosition(m4_position)

    def temperature(self):
        ''' Function for reading PT sensors 
        '''
        temp_vector = self._pt.getTemperature()
        return temp_vector
