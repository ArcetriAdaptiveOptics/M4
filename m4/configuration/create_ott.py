"""
Author(s)
---------
- Chiara Selmi: written in 2020
- Pietro Ferraiuolo: added DP in 2024
"""

class OTT:

    def __init__(
        self,
        parabola_slider,
        reference_mirror_slider,
        angle_rotator,
        parabola,
        reference_mirror,
        m4,
        dp,
        temperature_sensor,
        accelerometers,
    ):
        """The constructor"""
        self._parabola_slider = parabola_slider
        self._reference_mirror_slider = reference_mirror_slider
        self._angle_rotator = angle_rotator
        self._parabola = parabola
        self._reference_mirror = reference_mirror
        self._dp = dp if dp is not None else "DP not available"
        self._m4Exapode = m4
        self._pt = temperature_sensor
        self._acc = accelerometers

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
    def dp(self):
        return self._dp

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


    def __repr__(self):
        """
        Returns a string representation of the object.
        """
        dp_repr = ''*5+self.dp.__class__.__name__ if self.dp != "DP not available" else "DP not available"
        repr=\
f"""                  Optical Test Tower
---------------------------------------------------------
  .parabolaSlider         |          {self.parabolaSlider.__class__.__name__}
  .referenceMirrorSlider  |   {self.referenceMirrorSlider.__class__.__name__}
  .angleRotator           |            {self.angleRotator.__class__.__name__}
  .parabola               |                {self.parabola.__class__.__name__}
  .referenceMirror        |         {self.referenceMirror.__class__.__name__}
  .dp                     |            {dp_repr}
  .m4Exapode              |               {self.m4Exapode.__class__.__name__}
  .temperature            |      {self.temperature.__class__.__name__}
  .accelerometers         |          {self.accelerometers.__class__.__name__}
"""
        return repr