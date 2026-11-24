'''
Authors
  - C. Selmi: written in 2020
  - P. Ferraiuolo: rewritten in 2026 for the placeholder simulator
'''
import unittest
from unittest.mock import patch

from m4 import simulator
from m4.devices.ott import create_ott

_DEVICES = (
    "dm", "interferometer", "accelerometers", "angleRotator", "m4Exapode",
    "m4NySwitch", "dp", "parSlider", "par", "rmSlider", "rm", "tempSensors",
)


def _all_simulated_config():
    return {"SYSTEM": {"simulated.devices": {d: True for d in _DEVICES}}}


@patch("m4.devices.ott._lyc", side_effect=_all_simulated_config)
class TestStart(unittest.TestCase):

    def testCreationWithSimulatedDevices(self, _):
        ott, dm, interf = create_ott()

        expected = {
            "parabolaSlider": "ParabolaSlider",
            "referenceMirrorSlider": "ReferenceMirrorSlider",
            "angleRotator": "AngleRotator",
            "parabola": "Parabola",
            "referenceMirror": "ReferenceMirror",
            "m4Exapode": "M4Exapode",
            "m4Switcher": "M4NyasmithSwitcher",
            "temperature": "TemperatureSensor",
            "accelerometers": "Accelerometer",
        }
        for attr, name in expected.items():
            device = getattr(ott, attr)
            self.assertIsInstance(device, simulator._BasePlaceholder, attr)
            self.assertEqual(str(device), name)
        self.assertEqual(ott.dp, "DP not available")
        self.assertIsInstance(dm, simulator._DmPlaceholder)
        self.assertIsInstance(interf, simulator._InterfPlaceholder)

    def testNoDmNoInterfReturnsOnlyOtt(self, _):
        ott = create_ott(no_dm=True, no_interf=True)
        self.assertIsInstance(ott.parabola, simulator._BasePlaceholder)

    def testNoInterf(self, _):
        ott, dm = create_ott(no_interf=True)
        self.assertIsInstance(dm, simulator._DmPlaceholder)


class TestPlaceholders(unittest.TestCase):

    def testBaseMethodsAreNoOps(self):
        par = simulator.Parabola()
        self.assertIsNone(par.getPosition())
        self.assertIsNone(par.setPosition([0, 0, 1, 0, 0, 0]))
        with self.assertRaises(AttributeError):
            par.get_shape

    def testDmPlaceholder(self):
        dm = simulator.M4AU()
        self.assertIsNone(dm.get_shape())
        self.assertIsNone(dm.set_shape([0.0]))
        with self.assertRaises(AttributeError):
            dm.getPosition

    def testInterfPlaceholder(self):
        interf = simulator.PhaseCam()
        self.assertIsNone(interf.acquire_map())
        with self.assertRaises(AttributeError):
            interf.get_shape

    def testUnknownDeviceRaises(self):
        with self.assertRaises(AttributeError):
            simulator.NotADevice
