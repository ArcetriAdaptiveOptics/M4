'''
Authors
  - C. Selmi: written in February 2021
'''
import unittest
from m4.utils.parabola_identification import ParabolIdent


class Test(unittest.TestCase):

#    @unittest.skip('Non so leggere il file')
    def testFiduciali(self):
        pi = ParabolIdent()
        image = pi._imaTest('./img_0000.fits')
        coef1, coef2 = pi.testZernikeOnPar(image)
