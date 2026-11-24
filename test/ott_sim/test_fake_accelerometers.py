'''
Authors
  - C. Selmi: written in 2020
'''

import os
import shutil
import tempfile
import unittest
from unittest.mock import patch

import h5py
from deprecated.m4.simulator.fake_accelerometers import FakeAccelerometers


class TestFakeAccelerometers(unittest.TestCase):

    def setUp(self):
        self.acc = FakeAccelerometers()
        self._tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._tmpdir, ignore_errors=True)

    @patch('deprecated.m4.simulator.fake_accelerometers.fold_name', unsafe=True)
    def testDataAcquisition(self, mock_fold_name):
        mock_fold_name.ACCELEROMETERS_ROOT_FOLDER = self._tmpdir

        tt = self.acc.acquireData(recording_seconds=1)
        fname = os.path.join(self._tmpdir, tt + '.h5')

        self.assertTrue(os.path.isfile(fname))
        with h5py.File(fname, 'r') as hf:
            self.assertIn('Accelerometers', hf)
            self.assertEqual(hf['Accelerometers'].shape[0],
                             hf.attrs['TIME'].shape[0])
