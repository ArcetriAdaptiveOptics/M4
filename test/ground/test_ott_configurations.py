'''
Authors
  - C. Selmi: written in 2022
'''
import unittest
import os
import mock
from test.helper_test_library import testDataRootDir

class TestCalc(unittest.TestCase):

    def setUp(self):
        self.ott = self.createOtt()

    def tearDown(self):
        del self.ott

    @mock.patch('m4.ground.read_data.readFits_data', autospec=True)
    @mock.patch('numpy.load', autospec=True)
    def createOtt(self,  mock_rd, mock_load):
        from m4.configuration.ott import create_ott
        ott, interf, dm = create_ott()#os.path.join(testDataRootDir(), 'base',
                                       #        'Configurations', 'testConf.yaml'))
        return ott


    @mock.patch('m4.ground.read_data.readFits_data', autospec=True)
    def testConfigurations(self, mock_rd):
        from m4.ground.ott_configurations import OttConfigurations
        oc = OttConfigurations(self.ott)

        oc.move_to_segment_view(1, True)
        segment_view, rm_in = oc.get_configuration()
        self.assertEqual(segment_view, True)
        self.assertEqual(rm_in, True)

        oc.move_to_segment_view(1, False)
        segment_view, rm_in = oc.get_configuration()
        self.assertEqual(rm_in, False)

        oc.move_to_central_view(True)
        segment_view, rm_in = oc.get_configuration()
        self.assertEqual(segment_view, False)
        self.assertEqual(rm_in, True)

        oc.move_to_central_view(False)
        segment_view, rm_in = oc.get_configuration()
        self.assertEqual(rm_in, False)
