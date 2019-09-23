'''
@author: cs
'''

import unittest



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