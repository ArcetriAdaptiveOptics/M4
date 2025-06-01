'''
Authors
  - C. Selmi: written in 2022
'''
import unittest
import os
import mock
import numpy as np
from test.helper_test_library import testDataRootDir

class TestCalc(unittest.TestCase):

    def setUp(self):
        self.dm = self.createDM()
        self.readConfigData()

    def tearDown(self):
        del self.dm

    @mock.patch('m4.ground.read_data.readFits_data', autospec=True)
    @mock.patch('numpy.load', autospec=True)
    def createDM(self,  mock_rd, mock_load):
        from m4.configuration.ott import create_ott
        ott, interf, dm = create_ott()#os.path.join(testDataRootDir(), 'base',
                                     #         'Configurations', 'testConf.yaml'))
        return dm

    def readConfigData(self):
        from m4.ground import read_data
        root_test = os.path.join(testDataRootDir(), 'base', 'ottSim',
                                              'MIRROR_System', '20170430')
        self.dm.m4pupil = read_data.readFits_data(os.path.join(root_test,
                                                        'm4_mech_pupil-bin2.fits'))
        self.dm.m4ima = self.dm.m4pupil * 0.
        self.dm.CapsensGain = np.load(os.path.join(root_test,
                                              'CapsensGain.npy'))
        self.dm.ifmat = read_data.readFits_data(os.path.join(root_test,
                                                   'if_sect4_rot-bin2.fits'))
        self.dm.ifidx = read_data.readFits_data(os.path.join(root_test,
                                                    'if_idx4_rot-bin2.fits'))
        self.dm.vmat = read_data.readFits_data(os.path.join(root_test, 'ff_v_matrix.fits'))

    def test_fake_dm(self):
        n_acts = self.dm.getNActs()
        comm2 = self.dm.get_shape()
        #comm = np.ones(self.dm.getNActs())
        #self.dm.setActsCommand(comm)
