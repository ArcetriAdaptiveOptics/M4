"""
Author(s)
---------
- Pietro Ferraiuolo: written in 2025
- GitHub Copilot: revised in 2026

Description
-----------
OTT adapter for the EELT M4 hexapod Knowledge System.

The integration point in M4 is intentionally small and mirrors
[`FakeM4Exapode`](m4/simulator/fake_m4_exapode.py): only position read/write in
millimetres, plus metre helpers used by the optical simulator. The control path
implemented here is derived from the console application shipped in
`EELT-Console`, but the GUI is not reused: M4 talks directly to the KS through
ZMQ multipart messages and protobuf-compatible binary payloads.
"""

from __future__ import annotations

import logging as _logging
import struct as _struct
import time as _time

import numpy as _np
import zmq as _zmq

_logger = _logging.getLogger(__name__)

_WIRE_VARINT = 0
_WIRE_64BIT = 1
_WIRE_LENGTH_DELIMITED = 2
_WIRE_32BIT = 5

_DEFAULT_FILTER_VERSION = 17

_TOPIC_HASHES = {
    "ks:hp:mirror_pos": 1904128233,
    "ks:hp:setpoint_dof": -248262858,
    "ks:hp:loop": -684421844,
    "ks:hp:stop_motion": 532458397,
    "ks:acknowledge_error": 129238901,
}

_TYPE_HASHES = {
    "ks:hp:mirror_pos": 849988993,
    "ks:hp:setpoint_dof": 1485496452,
    "ks:hp:loop": 1816247398,
    "ks:hp:stop_motion": 2092280555,
    "ks:acknowledge_error": -1555128981,
}


def _encode_varint(value: int) -> bytes:
    if value < 0:
        value &= (1 << 64) - 1

    encoded = bytearray()
    while True:
        next_byte = value & 0x7F
        value >>= 7
        if value:
            encoded.append(next_byte | 0x80)
        else:
            encoded.append(next_byte)
            return bytes(encoded)


def _decode_varint(buffer: bytes, offset: int) -> tuple[int, int]:
    value = 0
    shift = 0
    while True:
        if offset >= len(buffer):
            raise ValueError("Unexpected end of protobuf buffer")

        next_byte = buffer[offset]
        offset += 1
        value |= (next_byte & 0x7F) << shift

        if not (next_byte & 0x80):
            return value, offset

        shift += 7
        if shift > 63:
            raise ValueError("Invalid protobuf varint")


def _encode_key(field_number: int, wire_type: int) -> bytes:
    return _encode_varint((field_number << 3) | wire_type)


def _encode_sint64(value: int) -> bytes:
    zigzag = (value << 1) ^ (value >> 63)
    return _encode_varint(zigzag)


def _build_filter_frame(topic_hash: int, filter_version: int) -> bytes:
    return _struct.pack(">i", int(topic_hash)) + b"\x00\x00" + bytes([filter_version]) + b"\x00"


def _build_header_frame(
    *,
    topic_hash: int,
    type_hash: int,
    sample_seq: int,
    timestamp: float | None = None,
) -> bytes:
    now = _time.time() if timestamp is None else float(timestamp)
    header = bytearray()
    header.extend(_encode_key(1, _WIRE_64BIT))
    header.extend(_struct.pack("<d", now))
    header.extend(_encode_key(2, _WIRE_VARINT))
    header.extend(_encode_varint(int(topic_hash)))
    header.extend(_encode_key(3, _WIRE_VARINT))
    header.extend(_encode_varint(int(sample_seq)))
    header.extend(_encode_key(4, _WIRE_VARINT))
    header.extend(_encode_varint(int(type_hash)))
    header.extend(_encode_key(11, _WIRE_VARINT))
    header.extend(_encode_varint(0))
    header.extend(_encode_key(12, _WIRE_VARINT))
    header.extend(_encode_sint64(0))
    return bytes(header)


def _encode_packed_doubles(field_number: int, values) -> bytes:
    packed = _np.asarray(values, dtype=float)
    if packed.shape != (6,):
        raise ValueError("Hexapod payloads require six values")

    payload = _struct.pack("<6d", *packed.tolist())
    return _encode_key(field_number, _WIRE_LENGTH_DELIMITED) + _encode_varint(len(payload)) + payload


def _encode_bool(field_number: int, value: bool) -> bytes:
    return _encode_key(field_number, _WIRE_VARINT) + _encode_varint(int(bool(value)))


def _extract_packed_doubles(payload: bytes, field_number: int) -> _np.ndarray:
    offset = 0
    while offset < len(payload):
        key, offset = _decode_varint(payload, offset)
        current_field = key >> 3
        wire_type = key & 0x07

        if wire_type == _WIRE_VARINT:
            _, offset = _decode_varint(payload, offset)
            continue

        if wire_type == _WIRE_64BIT:
            offset += 8
            continue

        if wire_type == _WIRE_32BIT:
            offset += 4
            continue

        if wire_type != _WIRE_LENGTH_DELIMITED:
            raise ValueError(f"Unsupported protobuf wire type {wire_type}")

        size, offset = _decode_varint(payload, offset)
        data = payload[offset : offset + size]
        offset += size

        if current_field != field_number:
            continue

        if len(data) != 48:
            raise ValueError("Expected six packed doubles")

        return _np.frombuffer(data, dtype="<f8").copy()

    raise ValueError(f"Field {field_number} not found in payload")


def _as_position(values, *, name: str) -> _np.ndarray:
    position = _np.asarray(values, dtype=float)
    if position.shape != (6,):
        raise ValueError(f"{name} must contain six values")
    return position


class M4Hexapode:
    """OTT-facing adapter for the EELT M4 hexapod."""

    def __init__(
        self,
        sub_ip: str,
        sub_port: int,
        pub_port: int,
        *,
        filter_version: int = _DEFAULT_FILTER_VERSION,
    ) -> None:
        self.sub_ip = sub_ip
        self.sub_port = sub_port
        self.pub_port = pub_port
        self.filter_version = int(filter_version)

        self._position = _np.zeros(6)
        self._topic_sample_seq: dict[str, int] = {
            topic: 0 for topic in _TYPE_HASHES if topic != "ks:hp:mirror_pos"
        }

        self._mirror_pos_filter = _build_filter_frame(
            _TOPIC_HASHES["ks:hp:mirror_pos"], self.filter_version
        )

        self._context = _zmq.Context()
        self._sub = self._context.socket(_zmq.SUB)
        self._pub = self._context.socket(_zmq.PUB)

        self._sub.connect(f"tcp://{self.sub_ip}:{self.sub_port}")
        self._sub.setsockopt(_zmq.SUBSCRIBE, self._mirror_pos_filter)

        self._pub.bind(f"tcp://*:{self.pub_port}")

    def close(self) -> None:
        self._sub.close(linger=0)
        self._pub.close(linger=0)
        self._context.term()

    def getPosition(self) -> _np.ndarray:
        latest = None
        while True:
            try:
                latest = self._sub.recv_multipart(flags=_zmq.NOBLOCK)
            except _zmq.Again:
                break

        if latest is None:
            return self._position.copy()

        try:
            self._position = self._decode_mirror_pos(latest)
        except ValueError as exc:
            _logger.warning("Invalid mirror_pos sample: %s", exc)

        return self._position.copy()

    def setPosition(self, absolute_position) -> _np.ndarray:
        """
        Set the hexapod absolute actuators position. The command is composed as:
        (tx, ty, tz, rx, ry, rz), where Translation is in millimeters and 
        Rotation in arcseconds.
        
        Parameters
        ----------
        absolute_position: ArrayLike
            The target position for the hexapod actuators, specified as a list 
            or array of six values corresponding to (tx, ty, tz, rx, ry, rz). 
            Translation values should be in millimeters and rotation values in 
            arcseconds.
        
        Returns
        -------
        position: ArrayLike
            The current position of the hexapod actuators after sending the 
            command.
        """
        position_mm = _as_position(
            absolute_position, name="absolute_position"
        )
        payload = _encode_packed_doubles(1, position_mm * 1e-3)
        self._send_command("ks:hp:setpoint_dof", payload)
        return self.getPosition()

    def getPositionInM(self) -> _np.ndarray:
        return self.getPosition() * 1e-3

    def setPositionInM(self, absolute_position_in_m) -> _np.ndarray:
        return self.setPosition(_as_position(absolute_position_in_m, name="absolute_position_in_m") * 1e3) * 1e-3

    def setLoopEnabled(self, enabled: bool) -> None:
        self._send_command("ks:hp:loop", _encode_bool(1, enabled))

    def stop(self) -> None:
        self._send_command("ks:hp:stop_motion", _encode_bool(1, False))

    def acknowledgeError(self) -> None:
        self._send_command("ks:acknowledge_error", _encode_bool(1, False))

    def _send_command(self, topic: str, payload: bytes) -> None:
        topic_hash = _TOPIC_HASHES[topic]
        type_hash = _TYPE_HASHES[topic]
        sample_seq = self._next_sample_seq(topic)
        frames = [
            _build_filter_frame(topic_hash, self.filter_version),
            _build_header_frame(
                topic_hash=topic_hash,
                type_hash=type_hash,
                sample_seq=sample_seq,
            ),
            payload,
        ]
        self._pub.send_multipart(frames)

    def _next_sample_seq(self, topic: str) -> int:
        sample_seq = self._topic_sample_seq[topic] + 1
        self._topic_sample_seq[topic] = sample_seq
        return sample_seq

    def _decode_mirror_pos(self, frames: list[bytes]) -> _np.ndarray:
        if len(frames) != 3:
            raise ValueError(f"Expected 3 frames, got {len(frames)}")

        if frames[0] != self._mirror_pos_filter:
            raise ValueError("Unexpected telemetry topic")

        position_m = _extract_packed_doubles(frames[2], 1)
        return position_m * 1e3
