"""
Author(s)
---------
- Chiara Selmi: written in 2020
                modified in 2022
- Pietro Ferraiuolo: modified in 2024
                     modified in 2026

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

```python
import m4
ott, dm, interf = m4.create_ott()
```

"""

try:
    import playsound # type: ignore
except:
    print("Playsound not installed, I'll be quiet...")

# from .opc_ua_controller import OpcUaController
from opticalib.core.config import load_yaml_config as _lyc
from ..configuration.ott_parameters import (
    EeltExapodeParameters as _ExapodePar,
    Sound as _sound
)

_OPCUA_DEVICES = ("parSlider", "rmSlider", "angleRotator", "par", "rm", "tempSensors")

# ------------------------------------------------ #
#   Check if we are running on the M4 development  #
#      platform and set up logging accordingly     #
# ------------------------------------------------ #
import os
_wsname = os.uname()[1]
if _wsname == 'm4dp': #TODO: right now the MICWS is `fedora` - need a custom name
    import Microgate.utils.setupLog as setupLog # type: ignore
    setupLog.consoleProfile()
 

def create_ott(
    *, no_dm: bool = False, no_interf: bool = False
) -> tuple[object, object, object]:
    """
    This function creates and initialize the OTT, creating all the devices, fake
    or real, accordingly to what specified in the .yaml configuration file.

    Parameters
    ----------
    no_dm: bool, optional
        If True, the Deformable Mirror will not be created, by default False
    no_interf: bool, optional
        If True, the Interferometer will not be created, by default False

    Returns
    -------
    ott: object
        The Optical Test Tower, comprehensive of:
        - the Parabola actuators and slider
        - the Reference Mirror actuators and slider
        - the Angle Rotator
        - M4's Exapode (or DP motors) and Nyasmith Switcher
        - Temperature Sensors
        - Accelerometers
    dm: object
        The M4 Deformable Mirror (if available)
    interf: object
        The 4D PhaseCam6110 Interferometer (if available)
    """
    config = _lyc()["SYSTEM"]["simulated.devices"]

#   TODO: Check if only one instance works

    if any(config[d] is False for d in _OPCUA_DEVICES):
        from .opc_ua_controller import OpcUaController

        opcUa = OpcUaController()

    # ------------------ #
    #       Sliders      #
    # ------------------ #

    # Parabola Slider
    if config["parSlider"] is True:
        from ..simulator import ParabolaSlider
        parabola_slider = ParabolaSlider()
    else:
        from .parabola_slider import OpcUaParabolaSlider

        #opcUa = OpcUaController()
        parabola_slider = OpcUaParabolaSlider(opcUa)

    # Reference Mirror Slider
    if config["rmSlider"] is True:
        from ..simulator import ReferenceMirrorSlider
        reference_mirror_slider = ReferenceMirrorSlider()
    else:
        from .reference_mirror_slider import OpcUaReferenceMirrorSlider

        #opcUa = OpcUaController()
        reference_mirror_slider = OpcUaReferenceMirrorSlider(opcUa)

    # Rotator
    if config["angleRotator"] is True:
        from ..simulator import AngleRotator
        angle_rotator = AngleRotator()
    else:
        from .angle_rotator import OpcUaAngleRotator

        #opcUa = OpcUaController()
        angle_rotator = OpcUaAngleRotator(opcUa)

    # ---------------------- #
    #   Alignment Controls   #
    # ---------------------- #

    # Parabola
    if config["par"] is True:
        from ..simulator import Parabola
        parab = Parabola()
    else:
        from .parabola import OpcUaParabola

        #opcUa = OpcUaController()
        parab = OpcUaParabola(opcUa)

    # Reference Mirror
    if config["rm"] is True:
        from ..simulator import ReferenceMirror
        reference_mirror = ReferenceMirror()
    else:
        from .reference_mirror import OpcUaReferenceMirror

        #opcUa = OpcUaController()
        reference_mirror = OpcUaReferenceMirror(opcUa)

    # DP
    if config["dp"] is True:
        dp = None
    else:
        from .dp_motors import ZmqDpMotors

        dp = ZmqDpMotors()

    # M4
    if config["m4Exapode"] is True:
        from ..simulator import M4Exapode
        m4hp = M4Exapode()
    else:
        try:
            from .eelt_exapode import M4Hexapode

            m4hp = M4Hexapode(
                sub_ip=_ExapodePar.sub_ip,
                sub_port=_ExapodePar.sub_port,
                pub_port=_ExapodePar.pub_port,
                filter_version=_ExapodePar.filter_version,
            )
        except Exception as e:
            print(f"M4 Hexapode unavailable: '{e}'")
            m4hp = None

    if config["m4NySwitch"] is True:
        from ..simulator import M4NyasmithSwitcher

        m4sm = M4NyasmithSwitcher()
    else:
        from .eelt_exapode import NyasmithSwitcher

        if config['m4Exapode'] is False and m4hp is not None:
            m4sm = NyasmithSwitcher.fromHexapode(m4hp)
        else:
            m4sm = NyasmithSwitcher(
                sub_ip=_ExapodePar.sub_ip,
                sub_port=_ExapodePar.sub_port,
                pub_port=_ExapodePar.pub_port,
                filter_version=_ExapodePar.filter_version,
            )


    # ------------------ #
    #       Others       #
    # ------------------ #

    # Temperature Sensors and Accelerometers
    if config["accelerometers"] is True:
        from ..simulator import Accelerometer
        accelerometers = Accelerometer()
    else:
        from .accelerometers import ZmqAccelerometers

        accelerometers = ZmqAccelerometers()

    if config["tempSensors"] is True:
        from ..simulator import TemperatureSensor
        temperature_sensor = TemperatureSensor()
    else:
        from .temperature_sensors import OpcUaTemperatureSensors

        #opcUa = OpcUaController()
        temperature_sensor = OpcUaTemperatureSensors(opcUa)

    # OTT creation
    ott = OTT(
        parabola_slider,
        reference_mirror_slider,
        angle_rotator,
        parab,
        reference_mirror,
        m4hp,
        m4sm,
        dp,
        temperature_sensor,
        accelerometers,
    )

    if _sound.PLAY is True:
        playsound.playsound(os.path.join(_sound.AUDIO_FILE_PATH, "ott-ini.mp3"))
        playsound.playsound(os.path.join(_sound.AUDIO_FILE_PATH, "ott-conf.mp3"))

    out = [ott]

    # ------------------ #
    #        DM          #
    # ------------------ #
    if not no_dm:
        try:
            if config["dm"] is True:
                from ..simulator import M4AU
                dm = M4AU()
            else:
                from opticalib.devices import DP
    
                dm = DP()
        except Exception as e:
            print(f"'{e}'")
            dm = None

        out.append(dm)

    # ------------------ #
    #   Interferometer   #
    # ------------------ #
    if not no_interf:
        if config["interferometer"] is True:
            from ..simulator import PhaseCam
            interf = PhaseCam()
        else:
            from opticalib import PhaseCam
            
            if _wsname=='m4dp':
                model = '6110_mic'
            else:
                model = '6110'
            interf = PhaseCam(model)

        out.append(interf)

    par0 = ott.parabolaSlider.getPosition()
    rm0  = ott.referenceMirrorSlider.getPosition()
    if par0 == 0:
        print('Sliding Tower not homed!')
    if rm0 == 0:
        print('Reference Mirror carriage not homed!')
    
    if len(out)==1:
        return out[0]

    return out


class OTT:

    def __init__(
        self,
        parabola_slider: object,
        reference_mirror_slider: object,
        angle_rotator: object,
        parabola: object,
        reference_mirror: object,
        m4hp: object,
        m4sm: object,
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
        self._m4Exapode = m4hp
        self._m4NySwitch = m4sm
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

    @property
    def m4Switcher(self):
        return self._m4NySwitch

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
