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
        self.client = mock.MagicMock()
        self.opc._client = self.client

    def tearDown(self):
        del self.opc

    @unittest.skip('Se non testi niente è inutile')
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

    def testStop(self):
        self.opc.stop()
        self.client.connect.assert_called_with()
        self.client.get_node.assert_called_with(self.opc.STOP_NODE)
        self.client.disconnect.assert_called_with()

    def testGetTemperatureVector(self):
        self.opc.get_temperature_vector()
        self.client.connect.assert_called_with()
        self.client.get_node.assert_called_with(self.opc.TEMPERATURE_NODE)
        self.client.disconnect.assert_called_with()

    def testGetPosition(self):
        self.opc.get_position(42)
        self.client.connect.assert_called_with()
        self.client.get_node.assert_called_with(
            "ns=7;s=MAIN.Drivers_input.f_PosAct[42]")
        self.client.disconnect.assert_called_with()

