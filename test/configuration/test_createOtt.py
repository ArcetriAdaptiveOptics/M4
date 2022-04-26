'''
Authors
  - C. Selmi: written in 2020
'''

import unittest
import numpy as np
from m4.ott_sim.fake_parabola_slider import FakeParabolaSlider
from m4.ott_sim.fake_reference_mirror_slider import FakeReferenceMirrorSlider
from m4.ott_sim.fake_angle_rotator import FakeAngleRotator
from m4.ott_sim.fake_parabola import FakeParabola
from m4.ott_sim.fake_reference_mirror import FakeReferenceMirror
from m4.ott_sim.fake_m4 import FakeM4
from m4.ott_sim.fake_temperature_sensors import FakeTemperatureSensors
from m4.configuration.create_ott import OTT
from numpy import testing
from m4.configuration.ott_parameters import OpcUaParameters
from m4.ott_sim.fake_accelerometers import FakeAccelerometers


class TestOtt(unittest.TestCase):

    def setUp(self):
        parabola_slider = FakeParabolaSlider()
        reference_mirror_slider = FakeReferenceMirrorSlider()
        angle_rotator = FakeAngleRotator()
        parabola = FakeParabola()
        reference_mirror = FakeReferenceMirror()
        m4 = FakeM4()
        temperature_sensor = FakeTemperatureSensors()
        accelerometers = FakeAccelerometers()
        self.ott = OTT(
            parabola_slider,
            reference_mirror_slider,
            angle_rotator,
            parabola,
            reference_mirror,
            m4,
            temperature_sensor,
            accelerometers)

    def testSlide(self):
        self.ott.parabolaSlider.setPosition(100)
        self.ott.referenceMirrorSlider.setPosition(-82.4)
        self.ott.angleRotator.setPosition(31.4)
        self.ott.parabola.setPosition(np.array([0, 0, -2, 3.3, -4, 0]))
        self.ott.referenceMirror.setPosition(np.array([0, 0, 0, 3.1, 4, 0]))
        self.ott.m4Exapode.setPosition(np.array([0, 0, 0, 33.3, 44.4, 0]))

        self.assertAlmostEqual(100,
                               self.ott.parabolaSlider.getPosition())
        self.assertAlmostEqual(-82.4,
                               self.ott.referenceMirrorSlider.getPosition())
        self.assertAlmostEqual(31.4,
                               self.ott.angleRotator.getPosition())
        testing.assert_allclose(self.ott.parabola.getPosition(),
                                np.array([0, 0, -2, 3.3, -4, 0]))
        testing.assert_allclose(self.ott.referenceMirror.getPosition(),
                                np.array([0, 0, 0, 3.1, 4, 0]))
        testing.assert_allclose(self.ott.m4Exapode.getPosition(),
                                np.array([0, 0, 0, 33.3, 44.4, 0]))
        testing.assert_allclose(self.ott.temperature.getTemperature(),
                                np.zeros(OpcUaParameters.num_PT_sensor))

