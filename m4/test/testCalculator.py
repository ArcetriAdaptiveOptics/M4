'''
@author: cs
'''

import unittest
import numpy as np

class TestCalc(unittest.TestCase):
    '''
    Contiene funzioni di test che utilizzano i dati presenti in
    /Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova
    '''

    def test_iffunctionsmaker_1(self):
        """ Test delle funzioni di influenza globali
        con command history casuale """
        print("IFFunctionsMaker test: \
            segment, global IF, shuffleCommandHistory")
        from m4.utils import create_device
        device = create_device.myDevice("segment")

        cmd_matrix_tag = 'matTestIF.fits'
        mode_vect_tag = 'vectTestIF.fits'
        amp_tag = 'ampTestIF.fits'
        n_push_pull = 3

        from m4 import sandbox
        tt = sandbox.testIFF_shuffleMeasureCreator(device, cmd_matrix_tag,
                                                   mode_vect_tag, amp_tag,
                                                   n_push_pull)
        an, prod, cube = sandbox.testIFF_an(tt)
        amp, sp_wf = sandbox.testIFF_spiano(an)
        aa = sp_wf.std()

        if aa < 1e-6:
            result = True
        else:
            result = False
        self.assertEqual(result, True, "Ohno")


    def test_iffunctionsmaker_2(self):
        """ Test delle funzioni di influenza zonali
        con command history ordinata """
        print("IFFunctionsMaker test: segment, zonal If, tidy command history")
        from m4.utils import create_device
        device = create_device.myDevice("segment")

        cmd_matrix_tag = 'matTestIF.fits'
        mode_vect_tag = 'vectTestIF.fits'
        amp_tag = 'ampTestIF.fits'
        n_push_pull = 3

        from m4 import sandbox
        tt = sandbox.testIFF_tidyMeasureCreator(device, cmd_matrix_tag,
                                                mode_vect_tag, amp_tag,
                                                n_push_pull)
        an, prod, cube = sandbox.testIFF_an(tt)
        amp, sp_wf = sandbox.testIFF_spiano(an)
        aa = sp_wf.std()

        if aa < 1e-6:
            result = True
        else:
            result = False
        self.assertEqual(result, True, "Ohno")


    def test_calibration_and_alignement(self):
        """test basato sui dati presenti in
        /Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/MixingIntMat/20190930_162714'
        l'immagine da allineare e le posizioni degli elementi ottici sono quelle in
        /Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/Allineamento/20191001_081344 """
        command_amp_vector = np.ones(5)*5.0e-06
        command_amp_vector[0] = 2.5e-05
        n_push_pull = 3
        print("Calibration test for PAR + RM")
        from m4.utils.optical_calibration import opt_calibration
        cal = opt_calibration()
        tt = cal.measureCalibrationMatrix(0, command_amp_vector, n_push_pull)
        mat, rec = cal.analyzerCalibrationMeasurement(tt, 3)
        print("Alignment test")
        from m4.utils.optical_alignment import opt_alignment
        al = opt_alignment(tt)
        cmd = al.opt_align()
        par, rm = al._testAlignment_loadInfoPositionFromFileFits(1)

        cmd_par = np.zeros(6)
        cmd_par[2] = cmd[0]
        cmd_par[3] = cmd[1]
        cmd_par[4] = cmd[2]
        cmd_rm = np.zeros(6)
        cmd_rm[3] = cmd[3]
        cmd_rm[4] = cmd[4]

        par[len(par)-1] = 0
        aa = par + cmd_par
        bb = rm + cmd_rm

        if (aa.sum() < 1e-10 and bb.sum() < 1e-10):
            result = True
        else:
            result = False
        self.assertEqual(result, True, "Ohno")
