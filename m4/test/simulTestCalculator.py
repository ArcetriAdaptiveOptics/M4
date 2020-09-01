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
        from m4.configuration.create_ott import OTT
        ott = OTT()
        from m4.alignment import Alignment
        a = Alignment(ott)
        tt = a.ott_calibration(command_amp_vector, n_push_pull, 3)

        ott.refflat(np.array([ 0.e+00,  0.e+00,  0.e+00,  1.e-07, -4.e-07,  0.e+00]))
        ott.parab(np.array([0,0,9.9999997e-06, 1.0836526e-07, 2.2898718e-09,0]))

        print("Alignment test")
        par_cmd, rm_cmd = a.ott_alignment(tt)

        aa = ott.parab() #nuova posizione
        bb = ott.refflat()

        if (aa.sum() < 1e-5 and bb.sum() < 1e-4):
            result = True
        else:
            result = False
        self.assertEqual(result, True, "Ohno")

    def test_calibration_and_alignement_m4(self):
        command_amp_vector = np.array([5.0e-06, 5.0e-06])
        n_push_pull = 1
        print("Calibration test for M4")
        from m4.configuration.create_ott import OTT
        ott = OTT()
        from m4.alignment import Alignment
        a = Alignment(ott)
        tt, zcoma, sur = a.m4_calibration(command_amp_vector, n_push_pull, 5)

        ott.m4(np.array([0,0,0,9.9999997e-06,0,0]))

        print("Alignment test")
        cmd = a.m4_alignment(zcoma, tt)