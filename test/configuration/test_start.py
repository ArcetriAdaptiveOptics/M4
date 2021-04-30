'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
from m4.configuration.start import create_ott
from m4.ott_sim.fake_parabola_slider import FakeParabolaSlider
from m4.ott_sim.fake_interferometer import FakeInterferometer


class FakeConfig():

    def __init__(self):
        self.simulated = 1


aFakeConfig = FakeConfig()


class TestStart(unittest.TestCase):

    def testCreationWithSimulatedDevices(self):
        aFakeConfig.simulated = 1
        ott, interf = create_ott(config=aFakeConfig)
        self.assertIsInstance(ott.parabolaSlider, FakeParabolaSlider)
        self.assertIsInstance(interf, FakeInterferometer)
