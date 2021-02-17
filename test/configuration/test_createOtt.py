'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import mock
import numpy as np

class TestCalc(unittest.TestCase):

    @mock.patch('m4.ground.read_data.readFits_data', autospec=True)
    @mock.patch('opcua.Client', autospec=True)
    def setUp(self, mock_readData, mock_opc):
        from m4.configuration.create_ott import OTT
        self.ott = OTT()
        self.ott.m4pupil = np.zeros((1237, 1237))
        self.ott.mask = np.zeros((512, 512))
        self.ott.parmask = np.zeros((512, 512))

    def tearDown(self):
        del self.ott

    @unittest.skip('Mancano i file per costruire il simulatore')
    def testTower(self):
        self.ott.slide()
        self.ott.slide(100)
        self.ott.rslide()
        self.ott.rslide(100)
        self.ott.angle()
        self.ott.angle(30)
        self.ott.parab()
        self.ott.parab(np.array([0,0,2,3,4,0]))
        self.ott.refflat()
        self.ott.refflat(np.array([0,0,0,3,4,0]))
        self.ott.m4()
        self.ott.m4(np.array([0,0,0,3,4,0]))
        #self.ott.temperature()
