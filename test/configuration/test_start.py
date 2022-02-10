'''
Authors
  - C. Selmi: written in 2020
'''
import os
import unittest
import mock
from m4.configuration.start import create_ott
from m4.ott_sim.fake_parabola_slider import FakeParabolaSlider
from m4.ott_sim.fake_interferometer import FakeInterferometer
from test.test_helper import testDataRootDir



class TestStart(unittest.TestCase):

    def testCreationWithSimulatedDevices(self):
        #aFakeConfig.simulated = 1
        ott, interf = create_ott(os.path.join(testDataRootDir(), 'base', 'Configurations', 'testConf.yaml'))
        self.assertIsInstance(ott.parabolaSlider, FakeParabolaSlider)
        self.assertIsInstance(interf, FakeInterferometer)

    #@mock.patch('m4.devices.interferometer.I4d4020', autospec=True)
    def testCreationWhitFakePathInYaml(self):
        ott, interf = create_ott(os.path.join(testDataRootDir(), 'base', 'Configurations', 'testConf2.yaml'))
