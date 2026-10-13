from m4.configuration.root import folders

import os as _os
import tempfile as _tempfile

from ruamel.yaml import YAML as _YAML

_gyml = _YAML()

with open(folders.CONFIGURATION_FILE, "r") as f:
    config = _gyml.load(f)

try:
    devices = [
        "dm",
        "interferometer",
        "accelerometers",
        "angleRotator",
        "m4Exapode",
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
        _tmp_path = None
        try:
            with _tempfile.NamedTemporaryFile(
                dir=_os.path.dirname(_os.path.abspath(folders.CONFIGURATION_FILE)),
                delete=False,
                mode="w",
            ) as _tmp:
                _tmp_path = _tmp.name
                _gyml.dump(config, _tmp)
            # keep the original file permissions (NamedTemporaryFile is 0600)
            _os.chmod(_tmp_path, _os.stat(folders.CONFIGURATION_FILE).st_mode & 0o7777)
            _os.replace(_tmp_path, folders.CONFIGURATION_FILE)
        except BaseException:
            if _tmp_path is not None and _os.path.exists(_tmp_path):
                _os.remove(_tmp_path)
            raise
        print("Updated configuration file with OTT Devices Settings.")
finally:
    for _n in ("config", "_YAML", "_gyml", "f", "devices", "device", "_system", "_simdev", "_missing", "_n", "_os", "_tempfile", "_tmp", "_tmp_path"):
        globals().pop(_n, None)


from m4.configuration.ott import create_ott
from .alignment import OttAligner

__all__ = [
    "create_ott",
    "OttAligner",
]
