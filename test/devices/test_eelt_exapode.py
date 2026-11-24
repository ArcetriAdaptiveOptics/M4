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


def _sm_pos_payload(actual, setpoint, actual_deg, setpoint_deg, commanded_deg):
    import struct

    payload = bytearray()
    payload += exapode._encode_key(1, exapode._WIRE_VARINT) + exapode._encode_varint(actual)
    payload += exapode._encode_key(2, exapode._WIRE_VARINT) + exapode._encode_varint(setpoint)
    for field, value in ((3, actual_deg), (4, setpoint_deg), (5, commanded_deg)):
        payload += exapode._encode_key(field, exapode._WIRE_64BIT) + struct.pack("<d", value)
    return bytes(payload)


class TestNyasmithSwitcher(unittest.TestCase):
    def _make(self, context_mock):
        sub_socket = MagicMock()
        pub_socket = MagicMock()
        sub_socket.recv_multipart.side_effect = exapode._zmq.Again()
        context = MagicMock()
        context.socket.side_effect = [sub_socket, pub_socket]
        context_mock.return_value = context
        device = exapode.NyasmithSwitcher("192.168.125.132", 56000, 57000)
        return device, sub_socket, pub_socket, context

    @patch("m4.devices.eelt_exapode._time.time", return_value=1234.5)
    @patch("m4.devices.eelt_exapode._zmq.Context", unsafe=True)
    def test_set_position_sends_console_compatible_frames(self, context_mock, _):
        device, sub_socket, pub_socket, _ = self._make(context_mock)

        returned = device.setPosition("nas_b")

        sub_socket.connect.assert_called_once_with("tcp://192.168.125.132:56000")
        sub_socket.setsockopt.assert_called_once_with(
            exapode._zmq.SUBSCRIBE, device._sm_pos_filter
        )
        pub_socket.bind.assert_called_once_with("tcp://*:57000")

        frames = pub_socket.send_multipart.call_args.args[0]
        topic = "ks:sm:setpoint_dof"
        self.assertEqual(
            frames[0],
            exapode._build_filter_frame(
                exapode._TOPIC_HASHES[topic], device.filter_version
            ),
        )
        self.assertEqual(
            frames[1],
            exapode._build_header_frame(
                topic_hash=exapode._TOPIC_HASHES[topic],
                type_hash=exapode._TYPE_HASHES[topic],
                sample_seq=1,
                timestamp=1234.5,
            ),
        )
        self.assertEqual(frames[2], b"\x08\x02")
        self.assertIs(returned, exapode.NasmythPosition.UNKNOWN_POS)

    @patch("m4.devices.eelt_exapode._zmq.Context", unsafe=True)
    def test_set_position_rejects_invalid_targets(self, context_mock):
        device, _, pub_socket, _ = self._make(context_mock)
        for bad in (0, "UNKNOWN_POS", 4, -1, "foo", None):
            with self.assertRaises(ValueError):
                device.setPosition(bad)
        pub_socket.send_multipart.assert_not_called()

    @patch("m4.devices.eelt_exapode._zmq.Context", unsafe=True)
    def test_commands_payloads_and_sequence(self, context_mock):
        device, _, pub_socket, _ = self._make(context_mock)

        with device.loopEnabled():
            pass
        device.stop()
        device.acknowledgeError()

        sent = [c.args[0] for c in pub_socket.send_multipart.call_args_list]
        expected_topics = [
            "ks:sm:loop",
            "ks:sm:loop",
            "ks:sm:stop_motion",
            "ks:acknowledge_error",
        ]
        for frames, topic in zip(sent, expected_topics):
            self.assertEqual(
                frames[0],
                exapode._build_filter_frame(exapode._TOPIC_HASHES[topic], 0),
            )
        self.assertEqual([f[2] for f in sent], [b"\x08\x01", b"\x08\x00", b"\x08\x00", b"\x08\x00"])
        self.assertEqual(device._topic_sample_seq["ks:sm:loop"], 2)

    @patch("m4.devices.eelt_exapode._zmq.Context", unsafe=True)
    def test_get_status_decodes_latest_sm_pos(self, context_mock):
        device, sub_socket, _, _ = self._make(context_mock)
        old = [device._sm_pos_filter, b"h", _sm_pos_payload(3, 3, 1.0, 2.0, 3.0)]
        new = [device._sm_pos_filter, b"h", _sm_pos_payload(1, 2, 12.5, 90.0, 30.25)]
        sub_socket.recv_multipart.side_effect = [old, new, exapode._zmq.Again()]

        status = device.getStatus()

        self.assertIs(status["actual"], exapode.NasmythPosition.NAS_A)
        self.assertIs(status["setpoint"], exapode.NasmythPosition.NAS_B)
        self.assertEqual(
            (status["actual_deg"], status["setpoint_deg"], status["commanded_deg"]),
            (12.5, 90.0, 30.25),
        )

    @patch("m4.devices.eelt_exapode._zmq.Context", unsafe=True)
    def test_invalid_samples_keep_last_value(self, context_mock):
        device, sub_socket, _, _ = self._make(context_mock)
        good = [device._sm_pos_filter, b"h", _sm_pos_payload(2, 2, 5.0, 5.0, 5.0)]
        sub_socket.recv_multipart.side_effect = [good, exapode._zmq.Again()]
        self.assertIs(device.getPosition(), exapode.NasmythPosition.NAS_B)

        bad_samples = [
            [device._sm_pos_filter, b"h", _sm_pos_payload(9, 2, 0.0, 0.0, 0.0)],
            [device._sm_pos_filter, b"h", _sm_pos_payload(1, 2, 0.0, 0.0, 0.0)[:-3]],
            [b"wrong", b"h", _sm_pos_payload(1, 1, 0.0, 0.0, 0.0)],
            [device._sm_pos_filter, b"h"],
        ]
        for bad in bad_samples:
            sub_socket.recv_multipart.side_effect = [bad, exapode._zmq.Again()]
            self.assertIs(device.getPosition(), exapode.NasmythPosition.NAS_B)
        sub_socket.recv_multipart.side_effect = exapode._zmq.Again()
        self.assertEqual(device.getPositionDeg(), 5.0)

    @patch("m4.devices.eelt_exapode._zmq.Context", unsafe=True)
    def test_from_hexapode_shares_publisher(self, context_mock):
        hp_sub, hp_pub, sm_sub = MagicMock(), MagicMock(), MagicMock()
        context = MagicMock()
        context.socket.side_effect = [hp_sub, hp_pub, sm_sub]
        context_mock.return_value = context

        hexapode = exapode.M4Hexapode("10.0.0.1", 56000, 57000, filter_version=17)
        switcher = exapode.NyasmithSwitcher.fromHexapode(hexapode)

        self.assertIs(switcher._pub, hp_pub)
        self.assertIs(switcher._context, context)
        hp_pub.bind.assert_called_once_with("tcp://*:57000")
        sm_sub.connect.assert_called_once_with("tcp://10.0.0.1:56000")
        self.assertEqual(switcher.filter_version, 17)

        switcher.stop()
        self.assertEqual(
            hp_pub.send_multipart.call_args.args[0][0],
            exapode._build_filter_frame(exapode._TOPIC_HASHES["ks:sm:stop_motion"], 17),
        )

        switcher.close()
        sm_sub.close.assert_called_once_with(linger=0)
        hp_pub.close.assert_not_called()
        context.term.assert_not_called()
