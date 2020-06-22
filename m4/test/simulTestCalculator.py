'''
Autors
  - C. Selmi: written in 2020
'''
import unittest
import numpy as np

class TestCalc(unittest.TestCase):

    def test_calibration_and_alignement_mix(self):
        command_amp_vector = np.array([5.0e-06, 5.0e-06, 2.5e-05, 5.0e-06, 5.0e-06])
        n_push_pull = 1
        print("Calibration test for PAR + RM")
        from m4.utils.optical_calibration import opt_calibration
        cal = opt_calibration()
        tt = cal.measureCalibrationMatrix(0, command_amp_vector, n_push_pull)
        mat, rec = cal.analyzerCalibrationMeasurement(tt, 2)

        from m4.configuration.create_ott import OTT
        ott = OTT()
        ott.slide(0.75)
        ott.angle(90.)
        ott.rslide(0.6)
        ott.refflat(np.array([ 0.e+00,  0.e+00,  0.e+00,  1.e-07, -4.e-07,  0.e+00]))
        ott.parab(np.array([0,0,9.9999997e-06, 1.0836526e-07, 2.2898718e-09,0]))

        print("Alignment test")
        from m4.utils.optical_alignment import opt_alignment
        al = opt_alignment(tt)
        par_cmd, rm_cmd = al.opt_align(ott)

        aa = ott.parab() + par_cmd
        bb = ott.refflat() + rm_cmd

        if (aa.sum() < 1e-10 and bb.sum() < 1e-10):
            result = True
        else:
            result = False
        self.assertEqual(result, True, "Ohno")
        