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
from test.helper_test_library import testDataRootDir



class TestStart(unittest.TestCase):

    @mock.patch('m4.ground.read_data.readFits_data', autospec=True)
    @mock.patch('numpy.load', autospec=True)
    def testCreationWithSimulatedDevices(self, mock_rd, mock_load):
        #aFakeConfig.simulated = 1
        ott, interf, dm = create_ott(os.path.join(testDataRootDir(), 'base', 'Configurations', 'testConf.yaml'))
        self.assertIsInstance(ott.parabolaSlider, FakeParabolaSlider)
        self.assertIsInstance(interf, FakeInterferometer)

    #@mock.patch('m4.devices.interferometer.I4d4020', autospec=True)
    @mock.patch('m4.ground.read_data.readFits_data', autospec=True)
    @mock.patch('numpy.load', autospec=True)
    def testCreationWhitFakePathInYaml(self, mock_rd, mock_load):
        ott, interf, dm = create_ott(os.path.join(testDataRootDir(), 'base', 'Configurations', 'testConf2.yaml'))
