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
        
        cmdMatrixTag= 'matTestIF.fits'
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
        
        cmdMatrixTag= 'matTestIF.fits'
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
        
        
    def testCalibrationAndAlignement(self):
        #test basato sui dati presenti in 
        #/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/MixingIntMat/20190930_162714'
        #l'immagine da allineare e le posizioni degli elementi ottici sono quelle in
        #/Users/rm/Desktop/Arcetri/M4/ProvaCodice/Immagini_prova/Allineamento/20191001_081344
        commandAmpVector=np.ones(5)*5.0e-06
        commandAmpVector[0]=2.5e-05 
        nPushPull=3
        print("Calibration test for PAR + RM")
        from m4.utils.opticalCalibration import Opt_Calibration
        cal= Opt_Calibration()
        tt= cal.measureCalibrationMatrix(0, commandAmpVector, nPushPull)
        mat, rec= cal.analyzerCalibrationMeasurement(tt, 3)
        print("Alignment test")
        from m4.utils.opticalAlignment import Opt_Alignment 
        al= Opt_Alignment(tt)
        cmd= al.opt_align()
        par, rm= al._testAlignment_loadInfoPositionFromFileFits(1)
        
        cmdPAR= np.zeros(6)
        cmdPAR[2]= cmd[0]
        cmdPAR[3]= cmd[1]
        cmdPAR[4]= cmd[2]
        cmdRM= np.zeros(6)
        cmdRM[3]= cmd[3]
        cmdRM[4]= cmd[4]
        
        par[len(par)-1]=0
        aa= par+cmdPAR
        bb= rm+cmdRM
        
        if (aa.sum()<1e-10 and bb.sum()<1e-10):
            result = True
        else:
            result = False
        self.assertEqual(result, True, "Ohno")

        

        
        
        
        
        
        
        
        