'''
Authors
  - C. Selmi: written in 2021
'''
import unittest
from types import SimpleNamespace
from unittest.mock import patch, MagicMock
import numpy as np
from m4.configuration.ott_parameters import OpcUaParameters


class TestOpcUaController(unittest.TestCase):
    """
    The opc ua client is mocked, so that no connection to the PLC is tried.
    The controller reads everything from the PLC telemetry node, whose content
    is faked here.
    """

    def setUp(self):
        from m4.devices.opc_ua_controller import OpcUaController
        self.telemetry = SimpleNamespace(
            RA=SimpleNamespace(Pos=SimpleNamespace(Act=42.5)),
            PAR=SimpleNamespace(
                PosAct=SimpleNamespace(PISTON=1.0, TIP=2.0, TILT=3.0)),
            System=SimpleNamespace(Temperatures=list(range(24))))
        self.nodes = {
            OpcUaParameters.NODE_ID_PLC_STATE: OpcUaParameters.PLC_RUN_STATE_VALUE,
            OpcUaParameters.NODE_ID_TELEMETRY: self.telemetry}

        def get_node(node_id):
            node = MagicMock()
            node.get_value.return_value = self.nodes[node_id]
            return node

        with patch('m4.devices.opc_ua_controller.Client') as client_cls, \
                patch('time.sleep'):  # skip the 0.5 s wait after first connect
            self.client = client_cls.return_value
            self.client.get_node.side_effect = get_node
            self.opc = OpcUaController()
        self.client.connect.reset_mock()
        self.client.disconnect.reset_mock()
        self.client.get_node.reset_mock()
        client_cls.assert_called_once_with(url=OpcUaParameters.server)

    def tearDown(self):
        del self.opc

    def testGetPosition(self):
        pos = self.opc.get_positions('RA')
        self.assertAlmostEqual(pos, 42.5)
        self.client.connect.assert_called_with()
        self.client.get_node.assert_any_call(OpcUaParameters.NODE_ID_TELEMETRY)
        self.client.disconnect.assert_called_with()

    def testGetPositionOfUnknownDeviceReturnsNone(self):
        self.assertIsNone(self.opc.get_positions('NOT_A_DEVICE'))

    def testGetActPositions(self):
        pos = self.opc.get_act_positions('PAR')
        self.assertEqual(pos, [1.0, 2.0, 3.0])
        self.client.disconnect.assert_called_with()

    def testGetTemperatureVector(self):
        temp = self.opc.get_temperature_vector()
        np.testing.assert_array_equal(temp, np.arange(24))
        self.client.connect.assert_called_with()
        self.client.get_node.assert_any_call(OpcUaParameters.NODE_ID_TELEMETRY)
        self.client.disconnect.assert_called_with()

    def testSetPositionRejectsTripodDevices(self):
        with self.assertRaises(ValueError):
            self.opc.set_position('PAR', 1.0)

    def testSetActPositionsValidatesInput(self):
        with self.assertRaises(ValueError):
            self.opc.set_act_positions('RA', [0, 0, 0])
        with self.assertRaises(ValueError):
            self.opc.set_act_positions('PAR', [0, 0])

    @unittest.skip("OpcUaController.stop() no longer exists: stop is issued "
                   "through send_command(mode='STOP', ...), which needs the "
                   "PLC command struct and is not testable without it")
    def testStop(self):
        self.opc.stop()
