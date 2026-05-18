import unittest
from unittest.mock import MagicMock, patch

import numpy as np
from numpy import testing

from m4.devices import eelt_exapode as exapode


class TestM4Exapode(unittest.TestCase):
    @patch("m4.devices.eelt_exapode._time.time", return_value=1234.5)
    @patch("m4.devices.eelt_exapode._zmq.Context", unsafe=True)
    def test_set_position_sends_console_compatible_frames(
        self, context_mock, _
    ):
        sub_socket = MagicMock()
        pub_socket = MagicMock()
        sub_socket.recv_multipart.side_effect = exapode._zmq.Again()

        context = MagicMock()
        context.socket.side_effect = [sub_socket, pub_socket]
        context_mock.return_value = context

        device = exapode.M4Exapode("192.168.125.132", 56000, 57000)
        returned = device.setPosition([0, 0, 1, 0, 0, 0])

        sub_socket.connect.assert_called_once_with("tcp://192.168.125.132:56000")
        sub_socket.setsockopt.assert_called_once_with(
            exapode._zmq.SUBSCRIBE, device._mirror_pos_filter
        )
        pub_socket.bind.assert_called_once_with("tcp://*:57000")

        sent_frames = pub_socket.send_multipart.call_args.args[0]
        self.assertEqual(
            sent_frames[0],
            exapode._build_filter_frame(
                exapode._TOPIC_HASHES["ks:hp:setpoint_dof"], 17
            ),
        )
        self.assertEqual(
            sent_frames[1],
            exapode._build_header_frame(
                topic_hash=exapode._TOPIC_HASHES["ks:hp:setpoint_dof"],
                type_hash=exapode._TYPE_HASHES["ks:hp:setpoint_dof"],
                sample_seq=1,
                timestamp=1234.5,
            ),
        )
        testing.assert_allclose(
            exapode._extract_packed_doubles(sent_frames[2], 1),
            np.array([0, 0, 0.001, 0, 0, 0]),
        )
        testing.assert_allclose(returned, np.zeros(6))

    @patch("m4.devices.eelt_exapode._zmq.Context", unsafe=True)
    def test_get_position_decodes_latest_actuals(self, context_mock):
        sub_socket = MagicMock()
        pub_socket = MagicMock()

        context = MagicMock()
        context.socket.side_effect = [sub_socket, pub_socket]
        context_mock.return_value = context

        device = exapode.M4Exapode("192.168.125.132", 56000, 57000)

        actual_m = np.array([0.001, -0.002, 0.003, 0.004, -0.005, 0.006])
        setpoint_m = np.zeros(6)
        frames = [
            device._mirror_pos_filter,
            b"header",
            exapode._encode_packed_doubles(1, actual_m)
            + exapode._encode_packed_doubles(2, setpoint_m),
        ]
        sub_socket.recv_multipart.side_effect = [frames, exapode._zmq.Again()]

        testing.assert_allclose(device.getPosition(), actual_m * 1e3)
