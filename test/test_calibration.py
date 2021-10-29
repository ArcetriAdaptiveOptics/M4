'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import numpy as np

class TestCalc(unittest.TestCase):

    def setUp(self):
        from m4.configuration import start
        self.ott, self._interf = start.create_ott()
        from m4.utils.optical_calibration import OpticalCalibration
        self.cal = OpticalCalibration(self.ott, self._interf)

    def tearDown(self):
        del(self.cal, self.ott)

    @unittest.skip('Salvataggio e lettura dati')
    def testCalibration(self):
        mixed_method = 1
        command_amp_vector = np.ones(5)
        n_push_pull = 1
        n_frames = 1
        mask_index = 4 #4 per il simulatore
        tt = self.cal.measureCalibrationMatrix(0, command_amp_vector,
                                               n_push_pull, n_frames, mixed_method)
        int_mat, rec = self.cal.analyzerCalibrationMeasurement(tt,
                                                                mask_index)
