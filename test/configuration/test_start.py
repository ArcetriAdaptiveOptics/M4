'''
Authors
  - C. Selmi: written in 2020
'''
import unittest
import mock
import numpy as np
from m4.configuration import config

class TestCalc(unittest.TestCase):

    @mock.patch('opcua.Client', autospec=True)
    @mock.patch('oaautils.i4d', autospect=True)
    def setUp(self, mock_opc, mock_i4d):
        from m4.configuration import start
        config.simulated = 1
        self.ott_sim, interf_sim = start.create_ott()
        config.simulated = 0
        self.ott, interf = start.create_ott()
#         self.ott.m4pupil = np.zeros((1237, 1237))
#         self.ott.mask = np.zeros((512, 512))
#         self.ott.parmask = np.zeros((512, 512))

    def tearDown(self):
        del self.ott, self.ott_sim

    def testTowerSimulator(self):
        self.ott_sim.slide()
        self.ott_sim.slide(100)
        self.ott_sim.rslide()
        self.ott_sim.rslide(100)
        self.ott_sim.angle()
        self.ott_sim.angle(30)
        self.ott_sim.parab()
        self.ott_sim.parab(np.array([0,0,2,3,4,0]))
        self.ott_sim.refflat()
        self.ott_sim.refflat(np.array([0,0,0,3,4,0]))
        self.ott_sim.m4()
        self.ott_sim.m4(np.array([0,0,0,3,4,0]))
        self.ott_sim.temperature()

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
        self.ott.temperature()

    def testInterfSimulator(self):
        pass

    def testInterfI4dArcetri(self):
        pass
