from astropy import config

from m4.configuration.root import folders

from opticalib.core.config import (
    load_yaml_config as _lyc,
    dump_yaml_config as _dyc,
)

config = _lyc(folders.CONFIGURATION_FILE)

try:
    devices = [
        "dm",
        "interferometer",
        "accelerometers",
        "angleRotator",
        "m4Exapode",
        "m4NySwitch",
        "dp",
        "parSlider",
        "par",
        "rmSlider",
        "rm",
        "tempSensors",
    ]
    _system = config["SYSTEM"]
    if "simulated.devices" not in _system or _system["simulated.devices"] is None:
        _system["simulated.devices"] = {}
    _simdev = _system["simulated.devices"]
    _missing = [device for device in devices if device not in _simdev]
    if _missing:
        for device in _missing:
            _simdev[device] = True
        _dyc(config, folders.CONFIGURATION_FILE)
        print("Updated configuration file with OTT Devices Settings.")
finally:
    for _n in (
        "config", "_lyc", "_dyc", "f", "devices", "device", "_system",
        "_simdev", "_missing", "_n",
    ):
        globals().pop(_n, None)


from m4.devices.ott import create_ott
from .alignment import OttAligner

__all__ = [
    "create_ott",
    "OttAligner",
]
