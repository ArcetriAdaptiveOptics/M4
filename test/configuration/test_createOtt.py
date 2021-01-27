'''
Authors
  - C. Selmi: written in 2020
'''
import unittest

class TestCalc(unittest.TestCase):

    def setUp(self):
        from m4.configuration.create_ott import OTT
        self.ott = OTT()

    def tearDown(self):
        del self.ott

    @unittest.skip('Mancano i file per costruire il simulatore')
    def testTower(self):
        self.ott.slide()
        self.ott.rslide()
        self.ott.angle()
        self.ott.parab()
        self.ott.refflat()
