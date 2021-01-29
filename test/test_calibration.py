'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import numpy as np

class TestCalc(unittest.TestCase):

    def setUp(self):
        from m4.configuration.create_ott import OTT
        self.ott = OTT()
        from m4.utils.optical_calibration import OpticalCalibration
        self.cal = OpticalCalibration()

    def tearDown(self):
        del(self.cal, self.ott)

    @unittest.skip('Salvataggio e lettura dati')
    def testCalibration(self):
        old_or_new = 1
        command_amp_vector = np.ones(5)
        n_push_pull = 1
        n_frames = 1
        mask_index = 4 #4 per il simulatore
        tt = self.cal.measureCalibrationMatrix(self.ott, 0, command_amp_vector,
                                                n_push_pull, n_frames, old_or_new)
        int_mat, rec = self.cal.analyzerCalibrationMeasurement(tt,
                                                                mask_index)
