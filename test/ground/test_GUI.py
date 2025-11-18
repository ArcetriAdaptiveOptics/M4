'''
Authors
  - C. Selmi: written in July 2021
'''
import os
import numpy as np
import unittest
from unittest.mock import patch
#from m4.gui import geometry_GUI

from m4.configuration import ott
from test.helper_test_library import testDataRootDir


class TestGui(unittest.TestCase):

    @unittest.skip('')
    @patch('m4.simulator.ott_images.conf', unsafe=True)
    @patch('m4.ground.read_data.readFits_data', unsafe=True)
    @patch('numpy.load', unsafe=True)
    def testGeometryGui(self, mock_conf, mock_rd, mock_load):
        want_mirror_root_folder = os.path.join(
            testDataRootDir(), 'base', 'ottSim',
            'MIRROR_System')
        mock_conf.MIRROR_FOLDER = want_mirror_root_folder

        ott, interf, dm = ott.create_ott()
        #g = GUI.Runner(ott)
        #.run()
