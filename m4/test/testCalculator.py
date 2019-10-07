'''
@author: cs
'''

import unittest
import numpy as np



class TestCalc(unittest.TestCase):

    def testIFFunctionsMaker1(self):
        print("IFFunctionsMaker test: segment, global IF, shuffle command history.")
        from m4.utils import createDevice 
        device= createDevice.myDevice("segment") 
        
        cmdMatrixTag= 'matShuffleTestIF.fits'
        modeVectTag= 'vectTestIF.fits'
        ampTag= 'ampTestIF.fits'
        nPushPull= 3
        
        from m4 import sandbox
        tt= sandbox.testIFF_shuffleMeasureCreator(device, cmdMatrixTag, modeVectTag, ampTag, nPushPull)
        an, prod, cube= sandbox.testIFF_an(tt)  
        amp, spWf= sandbox.testIFF_spiano(an)
        aa= spWf.std()
        
        if aa <1e-6:
            result = True
        else:
            result = False
        self.assertEqual(result, True, "Ohno")
        
      
    def testIFFunctionsMaker2(self):
        print("IFFunctionsMaker test: segment, zonal If, tidy command history")  
        from m4.utils import createDevice 
        device= createDevice.myDevice("segment")
        
        cmdMatrixTag= 'matTidyTestIF.fits'
        modeVectTag= 'vectTestIF.fits'
        ampTag= 'ampTestIF.fits'
        nPushPull= 3
        
        from m4 import sandbox
        tt= sandbox.testIFF_tidyMeasureCreator(device, cmdMatrixTag, modeVectTag, ampTag, nPushPull)
        an, prod, cube= sandbox.testIFF_an(tt)  
        amp, spWf= sandbox.testIFF_spiano(an)
        aa= spWf.std()
        
        if aa <1e-6:
            result = True
        else:
            result = False
        self.assertEqual(result, True, "Ohno")
        
        
    def testCalibration(self):
        #test basato sui dati presenti in 
        #/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/MixingIntMat/20190930_162714'
        commandAmpVector=np.ones(5)*5.0e-06
        nPushPull=3
        print("Calibration test for PAR + RM")
        from m4.opticalCalibration import Calibration
        cal= Calibration()
        tt= cal.measureCalibrationMatrix(0, commandAmpVector, nPushPull)
        mat, rec= cal.analyzerCalibrationMeasurement(tt)
        prod= np.dot(rec, mat) 
        np.fill_diagonal(prod,0) 
        if prod.sum()<1e-10:
            result = True
        else:
            result = False
        self.assertEqual(result, True, "Ohno")
        
        
        
        
        