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
        
        cmdMatrix=np.ones([892,800]) 
        cmdMatrix.T[30]=np.zeros(892)
        modeVect=np.array([10,11,12,13,14,15,16,17,18,19,30,31,32,33,34]) 
        amp= np.array([0,1,2,3,4,5,6,7,8,9,0.30,0.31,0.32,0.33,0.34])
    
        tt= IF.acq_IFFunctions(modeVect, 3, amp, cmdMatrix, 1)
        
        
        result = True
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
        modeVect=np.array([10,11,12,13,14,15,16,17,18,19,30,31,32,33,34]) 
        amp= np.array([0,1,2,3,4,5,6,7,8,9,0.30,0.31,0.32,0.33,0.34])
        
        tt= IF.acq_IFFunctions(modeVect, 3, amp, cmdMatrix)
        
        result = True
        self.assertEqual(result, True, "Ohno")