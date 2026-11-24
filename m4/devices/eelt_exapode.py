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

The module also provides `NyasmithSwitcher`, the adapter for the KS Nasmyth
switching mechanism (`ks:sm` topics), which selects one of the pre-set focal
stations in `NasmythPosition`.
"""

from __future__ import annotations

import logging as _logging
import struct as _struct
import time as _time
from contextlib import contextmanager as _ctm
from enum import IntEnum as _IntEnum

import numpy as _np
import zmq as _zmq

_logger = _logging.getLogger(__name__)

_WIRE_VARINT = 0
_WIRE_64BIT = 1
_WIRE_LENGTH_DELIMITED = 2
_WIRE_32BIT = 5

_DEFAULT_FILTER_VERSION = 0 # 17

_TOPIC_HASHES = {
    "ks:hp:mirror_pos": 1904128233,
    "ks:hp:setpoint_dof": -248262858,
    "ks:hp:loop": -684421844,
    "ks:hp:stop_motion": 532458397,
    "ks:acknowledge_error": 129238901,
    "ks:sm:pos": 1956015431,
    "ks:sm:setpoint_dof": 1195989849,
    "ks:sm:loop": 783282669,
    "ks:sm:stop_motion": 687566349,
}

_TYPE_HASHES = {
    "ks:hp:mirror_pos": 849988993,
    "ks:hp:setpoint_dof": 1485496452,
    "ks:hp:loop": 1816247398,
    "ks:hp:stop_motion": 2092280555,
    "ks:acknowledge_error": -1555128981,
    "ks:sm:pos": -1188135038,
    "ks:sm:setpoint_dof": 1527609140,
    "ks:sm:loop": 1462853566,
    "ks:sm:stop_motion": 1758558376,
}

_HP_COMMAND_TOPICS = (
    "ks:hp:setpoint_dof",
    "ks:hp:loop",
    "ks:hp:stop_motion",
    "ks:acknowledge_error",
)

_SM_COMMAND_TOPICS = (
    "ks:sm:setpoint_dof",
    "ks:sm:loop",
    "ks:sm:stop_motion",
    "ks:acknowledge_error",
)


class NasmythPosition(_IntEnum):
    """Admissible Nasmyth switcher positions (``NSPosition`` in ``m4ks.proto``)."""

    UNKNOWN_POS = 0
    NAS_A = 1
    NAS_B = 2
    PARK = 3


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


def _extract_scalar_fields(payload: bytes) -> dict[int, int | float | bytes]:
    """
    Decode the top-level fields of a protobuf payload.

    Varint fields are returned as ``int``, 64-bit fields as little-endian
    ``double``, 32-bit and length-delimited fields as raw ``bytes``. When a field
    number is repeated, the last occurrence wins (protobuf semantics for scalars).
    """
    fields: dict[int, int | float | bytes] = {}
    offset = 0
    while offset < len(payload):
        key, offset = _decode_varint(payload, offset)
        field_number = key >> 3
        wire_type = key & 0x07

        if wire_type == _WIRE_VARINT:
            value, offset = _decode_varint(payload, offset)
        elif wire_type == _WIRE_64BIT:
            if offset + 8 > len(payload):
                raise ValueError("Unexpected end of protobuf buffer")
            (value,) = _struct.unpack_from("<d", payload, offset)
            offset += 8
        elif wire_type == _WIRE_32BIT:
            if offset + 4 > len(payload):
                raise ValueError("Unexpected end of protobuf buffer")
            value = payload[offset : offset + 4]
            offset += 4
        elif wire_type == _WIRE_LENGTH_DELIMITED:
            size, offset = _decode_varint(payload, offset)
            if offset + size > len(payload):
                raise ValueError("Unexpected end of protobuf buffer")
            value = payload[offset : offset + size]
            offset += size
        else:
            raise ValueError(f"Unsupported protobuf wire type {wire_type}")

        fields[field_number] = value
    return fields


def _as_nasmyth_position(position) -> NasmythPosition:
    if isinstance(position, str):
        try:
            return NasmythPosition[position.strip().upper()]
        except KeyError:
            raise ValueError(f"Unknown Nasmyth position '{position}'") from None
    try:
        return NasmythPosition(int(position))
    except (TypeError, ValueError):
        raise ValueError(f"Unknown Nasmyth position '{position}'") from None


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
            topic: 0 for topic in _HP_COMMAND_TOPICS
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

        self._loop_enabled = False


    def getPosition(self) -> _np.ndarray:
        """
        Gets the hexapod absolute actuators position. The returned position array
        is composed as:
        
        ```python
            (tx, ty, tz, rx, ry, rz)
        ```
        
        where Translation (tx,ty,tz) is in millimeters and Rotation (rx,ry,rz)
        in arcseconds.

        Returns
        -------
        position: ArrayLike
            The current position of the hexapod actuators after sending the 
            command.
        """
        latest = None
        while True:
            try:
                latest = self._sub.recv_multipart(flags=_zmq.NOBLOCK)
            except _zmq.Again:
                break

        if latest is None:
            _logger.warning("No mirror_pos sample received. Reading last stored value.")
            return self._position.copy()

        try:
            self._position = self._decode_mirror_pos(latest)
        except ValueError as exc:
            _logger.warning("Invalid mirror_pos sample: %s", exc)

        return self._position.copy()

    def setPosition(self, absolute_position: _np.ndarray) -> _np.ndarray:
        """
        Set the hexapod absolute actuators position. The command is composed as:

        ```python
            (tx, ty, tz, rx, ry, rz)
        ```

        where Translation (tx,ty,tz) is in millimeters and Rotation (rx,ry,rz)
        in arcseconds.

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

    # def getPositionInM(self) -> _np.ndarray:
    #     return self.getPosition() * 1e-3

    # def setPositionInM(self, absolute_position_in_m) -> _np.ndarray:
    #     """Set position with iunput in meters in meters."""
    #     return self.setPosition(
    #         _as_position(
    #             absolute_position_in_m, name="absolute_position_in_m"
    #         ) * 1e3) * 1e-3

    @_ctm
    def loopEnabled(self):
        """
        Context manager for loop enabled state.
        
        It switches the loop ON while in the context and switches it back OFF on
        exit.
        """
        try:
            self.setLoopEnabled(True)
            yield
        except Exception as ex:
            raise ex
        finally:
            self.setLoopEnabled(False)

    def isLoopEnabled(self) -> bool:
        """Returns the current loop enabled state."""
        return self._loop_enabled

    def setLoopEnabled(self, enabled: bool) -> None:
        self._send_command("ks:hp:loop", _encode_bool(1, enabled))
        self._loop_enabled = enabled

    def acknowledgeError(self) -> None:
        """Acknowledges any hexapod error."""
        self._send_command("ks:acknowledge_error", _encode_bool(1, False))

    def stop(self) -> None:
        """Halts the hexapod motion."""
        self._send_command("ks:hp:stop_motion", _encode_bool(1, False))

    def close(self) -> None:
        """
        Close the sockets and context to both the subscriber and publisher.
        """
        self._sub.close(linger=0)
        self._pub.close(linger=0)
        self._context.term()

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


class NyasmithSwitcher:
    """
    OTT-facing adapter for the EELT M4 Nasmyth switching mechanism (``ks:sm``).

    The switcher is a single rotational axis that can only be commanded to one of
    the pre-set focal stations listed in :class:`NasmythPosition`. Commands and
    telemetry travel through the same KS ZMQ channel used by
    :class:`M4Hexapode`.

    Note
    ----
    The KS subscribes to a single console publisher endpoint, so an
    ``M4Hexapode`` and a ``NyasmithSwitcher`` living in the same process cannot
    both bind ``pub_port``. In that case build the switcher with
    :meth:`fromHexapode`, which reuses the hexapod publisher.
    """

    def __init__(
        self,
        sub_ip: str,
        sub_port: int,
        pub_port: int,
        *,
        filter_version: int = _DEFAULT_FILTER_VERSION,
        _context: _zmq.Context | None = None,
        _pub: _zmq.Socket | None = None,
    ) -> None:
        self.sub_ip = sub_ip
        self.sub_port = sub_port
        self.pub_port = pub_port
        self.filter_version = int(filter_version)

        self._position = NasmythPosition.UNKNOWN_POS
        self._setpoint = NasmythPosition.UNKNOWN_POS
        self._position_deg = _np.nan
        self._setpoint_deg = _np.nan
        self._commanded_deg = _np.nan
        self._topic_sample_seq: dict[str, int] = {
            topic: 0 for topic in _SM_COMMAND_TOPICS
        }

        self._sm_pos_filter = _build_filter_frame(
            _TOPIC_HASHES["ks:sm:pos"], self.filter_version
        )

        self._owns_pub = _pub is None
        self._owns_context = _context is None
        self._context = _zmq.Context() if _context is None else _context
        self._sub = self._context.socket(_zmq.SUB)
        self._sub.connect(f"tcp://{self.sub_ip}:{self.sub_port}")
        self._sub.setsockopt(_zmq.SUBSCRIBE, self._sm_pos_filter)

        if _pub is None:
            self._pub = self._context.socket(_zmq.PUB)
            self._pub.bind(f"tcp://*:{self.pub_port}")
        else:
            self._pub = _pub

        self._loop_enabled = False

    @classmethod
    def fromHexapode(cls, hexapode: M4Hexapode) -> "NyasmithSwitcher":
        """
        Build a switcher sharing the ZMQ context and publisher of ``hexapode``.

        A dedicated subscriber is still opened towards the same KS endpoint.
        Closing the switcher leaves the shared publisher and context untouched:
        close the hexapode afterwards to release them.
        """
        return cls(
            hexapode.sub_ip,
            hexapode.sub_port,
            hexapode.pub_port,
            filter_version=hexapode.filter_version,
            _context=hexapode._context,
            _pub=hexapode._pub,
        )

    def getPosition(self) -> NasmythPosition:
        """
        Gets the actual position of the switching mechanism.

        Returns
        -------
        position: NasmythPosition
            One of the pre-set Nasmyth positions. If no telemetry has been
            received yet, ``NasmythPosition.UNKNOWN_POS`` is returned.
        """
        self._update()
        return self._position

    def getPositionDeg(self) -> float:
        """
        Gets the actual angular position of the switching mechanism.

        Returns
        -------
        position: float
            Actual position in degrees (``nan`` if no telemetry was received).
        """
        self._update()
        return self._position_deg

    def getStatus(self) -> dict:
        """
        Gets the full ``ks:sm:pos`` sample.

        Returns
        -------
        status: dict
            Keys: ``actual`` and ``setpoint`` (:class:`NasmythPosition`),
            ``actual_deg``, ``setpoint_deg`` and ``commanded_deg`` (degrees).
        """
        self._update()
        return {
            "actual": self._position,
            "setpoint": self._setpoint,
            "actual_deg": self._position_deg,
            "setpoint_deg": self._setpoint_deg,
            "commanded_deg": self._commanded_deg,
        }

    def setPosition(self, position: NasmythPosition | int | str) -> NasmythPosition:
        """
        Commands the switching mechanism to one of the pre-set positions.

        The command is non-blocking: the movement takes time, so the returned
        value is the last position reported by the KS, not the target.

        Parameters
        ----------
        position: NasmythPosition | int | str
            Target station: ``NAS_A`` (1), ``NAS_B`` (2) or ``PARK`` (3).
            Strings are matched case-insensitively against the enum names.

        Returns
        -------
        position: NasmythPosition
            The current position of the switching mechanism after sending the
            command.
        """
        target = _as_nasmyth_position(position)
        if target is NasmythPosition.UNKNOWN_POS:
            raise ValueError("UNKNOWN_POS is not a valid setpoint")
        self._send_command(
            "ks:sm:setpoint_dof",
            _encode_key(1, _WIRE_VARINT) + _encode_varint(int(target)),
        )
        return self.getPosition()

    @_ctm
    def loopEnabled(self):
        """
        Context manager for loop enabled state.

        It switches the loop ON while in the context and switches it back OFF on
        exit.
        """
        try:
            self.setLoopEnabled(True)
            yield
        finally:
            self.setLoopEnabled(False)

    def isLoopEnabled(self) -> bool:
        """Returns the current loop enabled state."""
        return self._loop_enabled

    def setLoopEnabled(self, enabled: bool) -> None:
        """Closes (``True``) or opens (``False``) the switching mechanism loop."""
        self._send_command("ks:sm:loop", _encode_bool(1, enabled))
        self._loop_enabled = enabled

    def acknowledgeError(self) -> None:
        """
        Acknowledges any KS error.

        The console's ``ks:sm:acknowledge_error`` topic has no registered hash,
        so the global ``ks:acknowledge_error`` topic is used, as for the
        hexapod.
        """
        self._send_command("ks:acknowledge_error", _encode_bool(1, False))

    def stop(self) -> None:
        """Halts the switching mechanism motion."""
        self._send_command("ks:sm:stop_motion", _encode_bool(1, False))

    def close(self) -> None:
        """
        Close the subscriber, and the publisher and context if not shared with
        an :class:`M4Hexapode`.
        """
        self._sub.close(linger=0)
        if self._owns_pub:
            self._pub.close(linger=0)
        if self._owns_context:
            self._context.term()

    def _update(self) -> None:
        latest = None
        while True:
            try:
                latest = self._sub.recv_multipart(flags=_zmq.NOBLOCK)
            except _zmq.Again:
                break

        if latest is None:
            _logger.warning("No sm_pos sample received. Reading last stored value.")
            return

        try:
            self._decode_sm_pos(latest)
        except ValueError as exc:
            _logger.warning("Invalid sm_pos sample: %s", exc)

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

    def _decode_sm_pos(self, frames: list[bytes]) -> None:
        if len(frames) != 3:
            raise ValueError(f"Expected 3 frames, got {len(frames)}")

        if frames[0] != self._sm_pos_filter:
            raise ValueError("Unexpected telemetry topic")

        fields = _extract_scalar_fields(frames[2])
        try:
            actual = NasmythPosition(fields.get(1, 0))
            setpoint = NasmythPosition(fields.get(2, 0))
        except ValueError:
            raise ValueError("Invalid NSPosition value") from None
        decoded_deg = []
        for field_number in (3, 4, 5):
            value = fields.get(field_number, 0.0)
            if not isinstance(value, float):
                raise ValueError(f"Field {field_number} is not a double")
            decoded_deg.append(value)

        self._position = actual
        self._setpoint = setpoint
        self._position_deg, self._setpoint_deg, self._commanded_deg = decoded_deg
