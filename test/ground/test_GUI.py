'''
Authors
  - C. Selmi: written in July 2021
'''
import os
import numpy as np
import unittest
import mock
from m4.gui import geometry_GUI
from m4.configuration import start
from test.test_helper import testDataRootDir


class TestGui(unittest.TestCase):

    @unittest.skip('')
    @mock.patch('m4.ott_sim.ott_images.conf', unsafe=True)
    @mock.patch('m4.ground.read_data.readFits_data', unsafe=True)
    @mock.patch('numpy.load', unsafe=True)
    def testGeometryGui(self, mock_conf, mock_rd, mock_load):
        want_mirror_root_folder = os.path.join(
            testDataRootDir(), 'base', 'ottSim',
            'MIRROR_System')
        mock_conf.MIRROR_FOLDER = want_mirror_root_folder

        conf = os.path.join(testDataRootDir(), 'base', 'Configurations','testConf.yaml')
        ott, interf, dm = start.create_ott(conf)
        #g = GUI.Runner(ott)
        #.run()
