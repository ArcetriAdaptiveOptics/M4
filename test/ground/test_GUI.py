'''
Authors
  - C. Selmi: written in July 2021
'''
import os
import unittest
import mock
from m4.ground import GUI
from m4.configuration import start
from test.test_helper import testDataRootDir


class TestGui(unittest.TestCase):

    #@unittest.skip('')
    @mock.patch('m4.ott_sim.ott_images.conf', autospect=True)
    def testMainGui(self, mock_conf):
        want_mirror_root_folder = os.path.join(
            testDataRootDir(), 'base', 'ottSim',
            'MIRROR_System')
        mock_conf.MIRROR_FOLDER = want_mirror_root_folder

        conf = os.path.join(testDataRootDir(), 'base', 'ottSim','testConf.yaml')
        ott, interf = start.create_ott(conf)
        GUI.main(ott, interf)
