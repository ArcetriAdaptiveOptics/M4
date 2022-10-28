'''
Authors
  - C. Selmi: written in February 2021
'''
import unittest
import os
from m4.utils.parabola_identification import ParabolaCirculaPupilFromZernike

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'img_0000.fits')

class Test(unittest.TestCase):

    @unittest.skip('Non funziona il linalg.eig (non va da github)')
    def testFiduciali(self):
        pz = ParabolaCirculaPupilFromZernike()
        image = pz._imaTest(TESTDATA_FILENAME)
        coef1, coef2 = pz.zernike_on_parabola(image)
