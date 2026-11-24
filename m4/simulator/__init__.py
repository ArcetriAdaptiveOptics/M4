"""
OTT PLACEHOLDER MODULE FOR SIMULATED DEVICES
============================================

This module serves as a placeholder for simulated devices in the M4 project.

It replaces the old M4-OTT simulator, which used to be a real OTT devices simulator.
Now it's used with legacy conventions to put placeholders in the OTT class when
some devices are not available or do not need to be used.
"""

__DEVICES_CLASSES__ = {
    "dm": "M4AU",
    "interferometer": "PhaseCam",
    "accelerometers": "Accelerometer",
    "angleRotator": "AngleRotator",
    "m4Exapode": "M4Exapode",
    "m4NySwitch": "M4NyasmithSwitcher",
    "dp": "DP",
    "parSlider": "ParabolaSlider",
    "par": "Parabola",
    "rmSlider": "ReferenceMirrorSlider",
    "rm": "ReferenceMirror",
    "tempSensors": "TemperatureSensor",
}

__BASE_PLH_FNC__ = [
    "getPosition",
    "setPosition",
    "getPositionInM",
    "setPositionInM",
    "acquireData",
    "getTemperature",
]

__BASE_DM_FNC__ = [
    "get_shape",
    "set_shape",
    "upload_cmd_history",
    "run_cmd_history",
]

__BASE_INTERF_FNC__ = [
    "acquire_map",
    "capture",
    "produce"
]

class _BasePlaceholder:
    def __init__(self, cls_name):
        self.name = cls_name

    def __getattr__(self, name):
        if name in __BASE_PLH_FNC__:
            print(f"Device is not real: {self.__repr__()}")
            return lambda *args, **kwargs: None
        raise AttributeError(f"{self.name} has no attribute {name}")

    def __repr__(self):
        return f"<Placeholder for {self.name}>"

    def __str__(self):
        return f"{self.name}"

class _DmPlaceholder(_BasePlaceholder):
    def __getattr__(self, name):
        if name in __BASE_DM_FNC__:
            print(f"Device is not real: {self.__repr__()}")
            return lambda *args, **kwargs: None
        raise AttributeError(f"{self.name} has no attribute {name}")

class _InterfPlaceholder(_BasePlaceholder):
    def __getattr__(self, name):
        if name in __BASE_INTERF_FNC__:
            print(f"Device is not real: {self.__repr__()}")
            return lambda *args, **kwargs: None
        raise AttributeError(f"{self.name} has no attribute {name}")

def __getattr__(cls_name):
    if not cls_name in __DEVICES_CLASSES__.values():
        raise AttributeError(f"module {__name__} has no attribute {cls_name}")

    match cls_name:
        case 'M4AU':
            placeholder_class = _DmPlaceholder
        case 'PhaseCam':
            placeholder_class = _InterfPlaceholder
        case _:
            placeholder_class = _BasePlaceholder
    return placeholder_class(cls_name)


__all__ = [str(__getattr__(x)) for x in __DEVICES_CLASSES__.values()]