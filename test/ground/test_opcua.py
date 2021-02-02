'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
import mock


class Test(unittest.TestCase):

    @mock.patch('opcua.Client', autospec=True)
    def setUp(self, mock_client):
        from m4.ground.opc_ua_controller import OpcUaController
        self.opc = OpcUaController()

    def tearDown(self):
        del self.opc

    @unittest.skip('Funziona ma è lungo')
    def testName(self):
        self.opc.stop()
        self.opc.get_temperature_vector()
        self.opc.get_variables_positions()
        self.opc.get_position(0)
        self.opc.set_target_position(0, 1)
        self.opc.move_object(0)
        self.opc.readActsPositions(9, 10, 11)
        self.opc.setActsPositions(9, 10, 11, 1, 2, 3)
        self.opc.wait_for_stop(1)
