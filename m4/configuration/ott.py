"""
Author(s)
---------
- Chiara Selmi: written in 2020
                modified in 2022
- Pietro Ferraiuolo: modified in 2024

Description
-----------
Module which creates and/or connects to the M4's Optical Test Tower devices,
which are:
    - Accelerometers
    - Angle Rotator
    - M4's Exapode
    - DP alignment motors
    - Parabola (with parabola slider)
    - Reference Mirror (with reference mirror slider)
    - Temperature Sensors
    - Interferometer
    - Deformable Mirror
    
How to Use it
-------------
    >>> from m4.configuration import start
    >>> ott = start.create_ott()
"""

import os
import playsound
from m4.devices.dp_motors import ZmqDpMotors
from m4.devices.parabola import OpcUaParabola
from m4.devices.m4_exapode import OpcUaM4Exapode
from m4.simulator.fake_parabola import FakeParabola
from m4.configuration.ott_parameters import Sound
from m4.simulator.fake_m4_exapode import FakeM4Exapode
from m4.devices.angle_rotator import OpcUaAngleRotator
from m4.devices.accelerometers import ZmqAccelerometers
from m4.devices.opc_ua_controller import OpcUaController
from m4.simulator.fake_angle_rotator import FakeAngleRotator
from m4.devices.parabola_slider import OpcUaParabolaSlider
from m4.devices.reference_mirror import OpcUaReferenceMirror
from m4.simulator.fake_accelerometers import FakeAccelerometers
from m4.simulator.fake_interferometer import FakeInterferometer
from m4.simulator.fake_parabola_slider import FakeParabolaSlider
from m4.simulator.fake_reference_mirror import FakeReferenceMirror
from m4.devices.temperature_sensors import OpcUaTemperatureSensors
from m4.simulator.fake_deformable_mirror import FakeDeformableMirror
from m4.simulator.fake_temperature_sensors import FakeTemperatureSensors
from m4.devices.reference_mirror_slider import OpcUaReferenceMirrorSlider
from m4.simulator.fake_reference_mirror_slider import FakeReferenceMirrorSlider
from opticalib.core.read_config import load_yaml_config


def create_ott() -> object:
    """
    This function creates and initialize the OTT, creating all the devices, fake
    or real, accordingly to what specified in the .yaml configuration file.

    Returns
    -------
    ott: object
        The Optical Test Tower, comprehensive of:
         - the Parabola actuators and slider
         - the Reference Mirror actuators and slider
         - the Angle Rotator
         - M4's Exapode (or DP motors)
    """
    config = load_yaml_config()["SYSTEM"]["simulated.devices"]

    ###########
    # Sliders #
    ###########

    # Parabola Slider
    if config["parSlider"] is True:
        parabola_slider = FakeParabolaSlider()
    else:
        opcUa = OpcUaController()
        parabola_slider = OpcUaParabolaSlider(opcUa)

    # Reference Mirror Slider
    if config["rmSlider"] is True:
        reference_mirror_slider = FakeReferenceMirrorSlider()
    else:
        opcUa = OpcUaController()
        reference_mirror_slider = OpcUaReferenceMirrorSlider(opcUa)

    # Rotator
    if config["angleRotator"] is True:
        angle_rotator = FakeAngleRotator()
    else:
        opcUa = OpcUaController()
        angle_rotator = OpcUaAngleRotator(opcUa)

    #############
    # alignment #
    #############

    # Parabola
    if config["par"] is True:
        parab = FakeParabola()
    else:
        opcUa = OpcUaController()
        parab = OpcUaParabola(opcUa)

    # Reference Mirror
    if config["rm"] is True:
        reference_mirror = FakeReferenceMirror()
    else:
        opcUa = OpcUaController()
        reference_mirror = OpcUaReferenceMirror(opcUa)

    # DP
    if config["dp"] is True:
        dp = None
    else:
        dp = ZmqDpMotors()

    # M4
    if config["m4Exapode"] is True:
        m4 = FakeM4Exapode()
    else:
        opcUa = OpcUaController()
        m4 = OpcUaM4Exapode(opcUa)

    ##########
    # Others #
    ##########

    # Temperature Sensors and Accelerometers
    if config["accelerometers"] is True:
        accelerometers = FakeAccelerometers()
    else:
        accelerometers = ZmqAccelerometers()

    if config["tempSensors"] is True:
        temperature_sensor = FakeTemperatureSensors()
    else:
        opcUa = OpcUaController()
        temperature_sensor = OpcUaTemperatureSensors(opcUa)

    # OTT creation
    ott = OTT(
        parabola_slider,
        reference_mirror_slider,
        angle_rotator,
        parab,
        reference_mirror,
        m4,
        dp,
        temperature_sensor,
        accelerometers,
    )

    if Sound.PLAY is True:
        playsound.playsound(os.path.join(Sound.AUDIO_FILE_PATH, "ott-ini.mp3"))
        playsound.playsound(os.path.join(Sound.AUDIO_FILE_PATH, "ott-conf.mp3"))
    
    if config['dm'] is True:
        dm = FakeDeformableMirror()
    else:
        from opticalib import AdOpticaDm
        dm = AdOpticaDm()
    
    if config['interferometer'] is True:
        interf = FakeInterferometer(ott, dm)
    else:
        from opticalib import PhaseCam
        interf = PhaseCam('6110')

    return ott, dm, interf


class OTT:

    def __init__(
        self,
        parabola_slider: object,
        reference_mirror_slider: object,
        angle_rotator: object,
        parabola: object,
        reference_mirror: object,
        m4: object,
        dp: object,
        temperature_sensor: object,
        accelerometers: object,
    ):
        """The constructor"""
        self._parabola_slider = parabola_slider
        self._reference_mirror_slider = reference_mirror_slider
        self._angle_rotator = angle_rotator
        self._parabola = parabola
        self._reference_mirror = reference_mirror
        self._dp = dp if dp is not None else "DP not available"
        self._m4Exapode = m4
        self._pt = temperature_sensor
        self._acc = accelerometers

    # Elements position
    @property
    def parabolaSlider(self):
        return self._parabola_slider

    @property
    def referenceMirrorSlider(self):
        return self._reference_mirror_slider

    @property
    def angleRotator(self):
        return self._angle_rotator

    # Elements alignment
    @property
    def parabola(self):
        return self._parabola

    @property
    def referenceMirror(self):
        return self._reference_mirror

    @property
    def dp(self):
        return self._dp

    @property
    def m4Exapode(self):
        return self._m4Exapode

    # Other
    @property
    def temperature(self):
        return self._pt

    @property
    def accelerometers(self):
        return self._acc

    def __repr__(self):
        """
        Returns a string representation of the object.
        """
        dp_repr = (
            "" * 5 + self.dp.__class__.__name__
            if self.dp != "DP not available"
            else "DP not available"
        )
        repr = f"""                  Optical Test Tower
---------------------------------------------------------
  .parabolaSlider         |          {self.parabolaSlider.__class__.__name__}
  .referenceMirrorSlider  |   {self.referenceMirrorSlider.__class__.__name__}
  .angleRotator           |            {self.angleRotator.__class__.__name__}
  .parabola               |                {self.parabola.__class__.__name__}
  .referenceMirror        |         {self.referenceMirror.__class__.__name__}
  .dp                     |            {dp_repr}
  .m4Exapode              |               {self.m4Exapode.__class__.__name__}
  .temperature            |      {self.temperature.__class__.__name__}
  .accelerometers         |          {self.accelerometers.__class__.__name__}
"""
        return repr
