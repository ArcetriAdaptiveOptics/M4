"""
Author(s)
---------
    - Chiara Selmi: written in 2020
"""

class OTT():

    def __init__(self, parabola_slider, reference_mirror_slider, angle_rotator,
                 parabola, reference_mirror, m4, temperature_sensor, accelerometers):
        """The constructor """
        self._parabola_slider = parabola_slider
        self._reference_mirror_slider = reference_mirror_slider
        self._angle_rotator = angle_rotator
        self._parabola = parabola
        self._reference_mirror = reference_mirror
        self._m4Exapode = m4
        self._pt = temperature_sensor
        self._acc = accelerometers

#     @property
#     def slide(self):
#         return self._parabola_slider.getPosition()
#     @slide.setter
#     def slide(self, par_slider_position):
#         aa = self._parabola_slider.setPosition(par_slider_position)
#         print(aa)

# Elements position
    @property
    def parabolaSlider(self):
        return self._parabola_slider

    @property
    def referenceMirrorSlider(self):
        return self._reference_mirror_slider

    @property
    def angleRotator(self):
        return self._angle_rotator


# Elements alignment
    @property
    def parabola(self):
        return self._parabola

    @property
    def referenceMirror(self):
        return self._reference_mirror

    @property
    def m4Exapode(self):
        return self._m4Exapode

# Other
    @property
    def temperature(self):
        return self._pt

    @property
    def accelerometers(self):
        return self._acc
