'''
Authors
  - C. Selmi: written in 2020
'''

import os
import unittest
from unittest.mock import patch
from m4.simulator.fake_accelerometers import FakeAccelerometers
from test.helper_test_library import testDataRootDir


class TestFakeAccelerometers(unittest.TestCase):

    def setUp(self):
        self.acc = FakeAccelerometers()

    def tearDown(self):
        os.unlink(self._fname)

    @patch('m4.simulator.fake_accelerometers.fold_name', unsafe=True)
    def testDataAcquisition(self, mock_fold_name):
        want_acc_root_folder = os.path.join(
            testDataRootDir(), 'base', 'M4Data',
            'OPTData', 'AccelerometersData')
        mock_fold_name.ACC_ROOT_FOLDER = want_acc_root_folder

        tt = self.acc.acquireData()
        name = tt + '.h5'
        self._fname = os.path.join(want_acc_root_folder, name)
        print("filename %s " % self._fname)
        # TODO: test that the fname contains the correct matrix and attrs
