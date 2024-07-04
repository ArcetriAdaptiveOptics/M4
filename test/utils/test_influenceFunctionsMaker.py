'''
Authors
  - C. Selmi:  written in 2022
'''
import unittest
import os
import numpy as np
from m4.utils.influence_functions_maker import IFFunctionsMaker
from test.helper_test_library import testDataRootDir
import mock

class TestInfluenceFunctionsMaker(unittest.TestCase):

    def setUp(self):
        self.dm = self._createDeformableMirror()
        self.interf = self._createInterferometer()

    def tearDown(self):
        self.dm
        self.interf

    @mock.patch('astropy.io.open', autospec=None)
    @mock.patch('m4.utils.influence_functions_maker.IFFunctionsMaker._storageFolder', autospec=True)
    @mock.patch('astropy.io.fits.writeto', autospec=None)
    @mock.patch('m4.type.modalAmplitude.ModalAmplitude._storageFolder', autospec=True)
    @mock.patch('m4.type.modalBase.ModalBase._storageFolder', autospec=True)
    @mock.patch('m4.ground.tracking_number_folder.os.makedirs', autospec=True)
    @mock.patch('m4.type.commandHistory.CmdHistory.saveInfo', autospec=True)
    def testIFFsAcquisition(self, mock, mockcmd, mockFilepath1, mockFilepath2, mock_fits,
                            mock_folder, mock_open):
        modalBaseTag = 'Hadamard3.fits'
        ampTag = 'ampTest3.fits'
        iff = IFFunctionsMaker(self.dm, self.interf)
        iff._interf.save_phasemap = self._skipSave
        iff._interf.acquire_phasemap = self._skipAcq
        n_push_pull = 1
        #inserire mock skip salvataggio dati
        mockFilepath2.return_value = os.path.join(testDataRootDir(), 'base', 'M4Data', 'OPTData', 'ModalAmplitude')
        mockFilepath1.return_value = os.path.join(testDataRootDir(), 'base', 'M4Data', 'OPTData', 'ModalBase')
        tt = iff.acq_IFFunctions(n_push_pull, ampTag, modalBaseTag,
                                 shuffle=False, template=None)
        tt2 = iff.acq_IFFunctions(n_push_pull, ampTag, modalBaseTag,
                                 shuffle=True, template=np.array([1,-1,1]))

        mock_folder.return_value = os.path.join(testDataRootDir(), 'base',
                                                  'M4Data/OPTData/IFFunctions')
        folder = os.path.join(testDataRootDir(), 'base','M4Data/OPTData/IFFunctions',
                              tt2)
        #iff._saveInfoFile(folder, 0)

    @mock.patch('m4.utils.influence_functions_maker.IFFunctionsMaker._storageFolder', autospec=True)
    def testReload(self, mock_folder):
        tt = '20220309_142454'
        mock_folder.return_value = os.path.join(testDataRootDir(), 'base',
                                                  'M4Data/OPTData/IFFunctions')
        IFFunctionsMaker.loadInfo(tt, 0)

    def _createDeformableMirror(self):
        dm = DMtest()
        return dm

    @mock.patch('m4.ground.read_data.readFits_data', autospec=True)
    @mock.patch('numpy.load', autospec=True)
    def _createInterferometer(self, mock_rd, mock_load):
        from m4.configuration.start import create_ott
        from m4.ott_sim.fake_interferometer import FakeInterferometer
        ott, interf, dm = create_ott()#os.path.join(testDataRootDir(), 'base', 'Configurations', 'testConf.yaml'))
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
        acts = 3
        return acts
    def setActsCommand(self, cmd):
        pass
