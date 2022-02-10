'''
Authors
  - C. Selmi:  written in 2022
'''
import unittest
import os
from m4.utils.influence_functions_maker import IFFunctionsMaker
from test.test_helper import testDataRootDir
import mock

class TestInfluenceFunctionsMaker(unittest.TestCase):

    def setUp(self):
        self.dm = self._createDeformableMirror()
        self.interf = self._createInterferometer()

    def tearDown(self):
        self.dm
        self.interf

    @mock.patch('m4.type.modalAmplitude.ModalAmplitude._storageFolder', autospec=True)
    @mock.patch('m4.type.modalBase.ModalBase._storageFolder', autospec=True)
    @mock.patch('m4.ground.tracking_number_folder.os.makedirs', autospec=True)
    @mock.patch('m4.type.commandHistory.CmdHistory.saveInfo', autospec=True)
    def testIFFsAcquisition(self, mock, mockcmd, mockFilepath1, mockFilepath2):
        modalBaseTag = 'Hadarmard10.fits'
        ampTag = 'ampTest10.fits'
        iff = IFFunctionsMaker(self.dm, self.interf)
        iff._interf.save_phasemap = self._skipSave
        iff._interf.acquire_phasemap = self._skipAcq
        n_push_pull = 1
        #inserire mock skip salvataggio dati
        mockFilepath2.return_value = os.path.join(testDataRootDir(), 'base', 'M4Data', 'OPTData', 'ModalAmplitude')
        mockFilepath1.return_value = os.path.join(testDataRootDir(), 'base', 'M4Data', 'OPTData', 'ModalBase')
        tt = iff.acq_IFFunctions(n_push_pull, ampTag, modalBaseTag,
                                 shuffle=False, template=None)

    def _createDeformableMirror(self):
        dm = DMtest()
        return dm

    #@mock.patch('m4.ott_sim.fake_interferometer.FakeInterferometer.save_phasemap', autospec=True)
    #@mock.patch('m4.ott_sim.fake_interferometer.FakeInterferometer.acquire_phasemap', autospec=True)
    def _createInterferometer(self):
        from m4.configuration.start import create_ott
        from m4.ott_sim.fake_interferometer import FakeInterferometer
        ott, interf = create_ott(os.path.join(testDataRootDir(), 'base', 'ottSim', 'testConf.yaml'))
        self.assertIsInstance(interf, FakeInterferometer)
        return interf

    def _skipAcq(self, aa):
        pass

    def _skipSave(self, a, b, c):
        pass


class DMtest():
    def __init__(self):
        pass
    def getNActs(self):
        acts = 10
        return acts
    def setActsCommand(self, cmd):
        pass
