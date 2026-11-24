"""
OTT PLACEHOLDER MODULE FOR SIMULATED DEVICES
============================================

This module serves as a placeholder for simulated devices in the M4 project.

It replaces the old M4-OTT simulator, which used to be a real OTT devices simulator.
Now it's used with legacy conventions to put placeholders in the OTT class when
some devices are not available or do not need to be used.
"""
_BASE_FNC = ("getPosition", "setPosition", "getPositionInM", "setPositionInM",
             "acquireData", "getTemperature")
_DM_FNC = ("get_shape", "set_shape", "upload_cmd_history", "run_cmd_history")
_INTERF_FNC = ("acquire_map", "capture", "produce")


class _BasePlaceholder:
    """Stateless stand-in for an unavailable device: whitelisted calls are no-ops."""

    _methods: frozenset[str] = frozenset(_BASE_FNC)

    def __getattr__(self, name):
        # only reached for attributes not found normally; dunders (pickle/copy
        # probes) must raise, never be faked
        if name.startswith("__") or name not in type(self)._methods:
            raise AttributeError(f"{type(self).__name__} placeholder has no attribute {name!r}")

        def _noop(*args, **kwargs):
            print(f"{type(self).__name__}.{name}() called on a placeholder: no-op")
            return None

        _noop.__name__ = name
        return _noop

    def __repr__(self):
        return f"<Placeholder for {type(self).__name__}>"

    def __str__(self):
        return type(self).__name__


class _DmPlaceholder(_BasePlaceholder):
    _methods = frozenset(_DM_FNC)


class _InterfPlaceholder(_BasePlaceholder):
    _methods = frozenset(_INTERF_FNC)


def _placeholder(name, base=_BasePlaceholder):
    return type(name, (base,), {"__module__": __name__, "__qualname__": name})


M4AU = _placeholder("M4AU", _DmPlaceholder)
PhaseCam = _placeholder("PhaseCam", _InterfPlaceholder)
Accelerometer = _placeholder("Accelerometer")
AngleRotator = _placeholder("AngleRotator")
M4Exapode = _placeholder("M4Exapode")
M4NyasmithSwitcher = _placeholder("M4NyasmithSwitcher")
DP = _placeholder("DP")
ParabolaSlider = _placeholder("ParabolaSlider")
Parabola = _placeholder("Parabola")
ReferenceMirrorSlider = _placeholder("ReferenceMirrorSlider")
ReferenceMirror = _placeholder("ReferenceMirror")
TemperatureSensor = _placeholder("TemperatureSensor")

__all__ = ["M4AU", "PhaseCam", "Accelerometer", "AngleRotator", "M4Exapode",
           "M4NyasmithSwitcher", "DP", "ParabolaSlider", "Parabola",
           "ReferenceMirrorSlider", "ReferenceMirror", "TemperatureSensor"]
