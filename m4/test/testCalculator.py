'''
@author: cs
'''

import unittest
import numpy as np 
import os
from m4.ground.configuration import Configuration


class TestCalc(unittest.TestCase):

    def testIFFunctionsMaker1(self):
        print("IFFunctionsMaker test: segment, global IF, shuffle command history.")
        from m4.utils import createDevice 
        device= createDevice.myDevice("segment") 
        from m4.influenceFunctionsMaker import IFFunctionsMaker
        IF= IFFunctionsMaker(device) 
        
        cmdMatrix=np.zeros([892,800]) 
        cmdMatrix.T[30]=np.ones(892)
        modeVect=np.arange(25)
        amp= np.arange(25)+1
        nPushPull= 3
        
        from m4 import sandbox
        tt= sandbox.testIFF_shuffleMeasureCreator(device, cmdMatrix, modeVect, amp, nPushPull)
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
        from m4.influenceFunctionsMaker import IFFunctionsMaker
        IF= IFFunctionsMaker(device)
        
        cmdMatrix=np.zeros((892,800))  
        for i in range(800):
            cmdMatrix[i][i]=1
        modeVect=np.arange(25)
        amp= np.arange(25)+1
        nPushPull= 3
        
        from m4 import sandbox
        tt= sandbox.testIFF_tidyMeasureCreator(device, cmdMatrix, modeVect, amp, nPushPull)
        an, prod, cube= sandbox.testIFF_an(tt)  
        amp, spWf= sandbox.testIFF_spiano(an)
        aa= spWf.std()
        
        if aa <1e-6:
            result = True
        else:
            result = False
        self.assertEqual(result, True, "Ohno")