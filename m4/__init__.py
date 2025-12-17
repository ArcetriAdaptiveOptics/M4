from m4.configuration.root import folders

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
    if not all([device in config["SYSTEM"]["simulated.devices"] for device in devices]):
        raise KeyError
except KeyError:
    config["SYSTEM"]["simulated.devices"]["dm"] = True
    config["SYSTEM"]["simulated.devices"]["interferometer"] = True
    config["SYSTEM"]["simulated.devices"]["accelerometers"] = True
    config["SYSTEM"]["simulated.devices"]["angleRotator"] = True
    config["SYSTEM"]["simulated.devices"]["m4Exapode"] = True
    config["SYSTEM"]["simulated.devices"]["dp"] = True
    config["SYSTEM"]["simulated.devices"]["parSlider"] = True
    config["SYSTEM"]["simulated.devices"]["par"] = True
    config["SYSTEM"]["simulated.devices"]["rmSlider"] = True
    config["SYSTEM"]["simulated.devices"]["rm"] = True
    config["SYSTEM"]["simulated.devices"]["tempSensors"] = True
    with open(folders.CONFIGURATION_FILE, "w") as f:
        _gyml.dump(config, f)
    print("Updated configuration file with OTT Devices Settings.")
finally:
    del config, _YAML, _gyml, f, devices


from m4.configuration.ott import create_ott
from .alignment import OttAligner
from .ottiff import OttIffAcquisition
from .spl_meas import SplAcquirer, SplAnalyzer


def create_spl():
    from opticalib.devices import AVTCamera
    from plico_motor import motor  # type: ignore

    # FIXME: would be better to read this from a configuration somewhere
    camera = AVTCamera(name="SplCam0")
    filter = motor("localhost", 7300, axis=1)

    return SplAcquirer(camera=camera, filter_motor=filter)


__all__ = [
    "create_ott",
    "create_spl",
    "OttAligner",
    "OttIffAcquisition",
    "SplAcquirer",
    "SplAnalyzer",
]
