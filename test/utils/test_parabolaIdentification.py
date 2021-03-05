'''
Authors
  - C. Selmi: written in February 2021
'''
import unittest
import os
from m4.utils.parabola_identification import ParabolIdent

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'img_0000.fits')

class Test(unittest.TestCase):


    def testFiduciali(self):
        pi = ParabolIdent()
        image = pi._imaTest(TESTDATA_FILENAME)
        coef1, coef2 = pi.testZernikeOnPar(image)
