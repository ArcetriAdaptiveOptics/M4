'''
Authors
  - C. Selmi: written in 2020
'''

import unittest
import mock
from m4.ott_sim.fake_accelerometers import FakeAccelerometers
from m4.configuration.ott_parameters import OpcUaParameters


class TestFakeAccelerometers(unittest.TestCase):

    def setUp(self):
        self.acc = FakeAccelerometers()

    @mock.patch('h5py.File', autospect=True)
    def testDataAcquisition(self, mock_h5py):
        self.acc.acquireData()
        #mock_h5py.attrs.assert_called_once_with(OpcUaParameters.accelerometers_dt)
        #mock_h5py.attrs.assert_called_once_with(OpcUaParameters.accelerometers_plc_id)
        #controllare le dimensioni di signal
