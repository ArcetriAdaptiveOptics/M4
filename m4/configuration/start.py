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
    >>> ott, interf, dm = start.create_ott(config_file_path)
"""

import os
import playsound
from m4.configuration.create_ott import OTT
from m4.devices.dp_motors import ZmqDpMotors
from m4.devices.interferometer import I4d6110
from m4.devices.parabola import OpcUaParabola
from m4.devices.m4_exapode import OpcUaM4Exapode
from m4.ott_sim.fake_parabola import FakeParabola
from m4.configuration.ott_parameters import Sound
from m4.ott_sim.fake_m4_exapode import FakeM4Exapode
from m4.devices.angle_rotator import OpcUaAngleRotator
from m4.ott_sim.fake_deformable_mirror import FakeM4DM
from m4.configuration import update_folder_paths as ufp
from m4.devices.accelerometers import ZmqAccelerometers
from m4.devices.opc_ua_controller import OpcUaController
from m4.ott_sim.fake_angle_rotator import FakeAngleRotator
from m4.devices.parabola_slider import OpcUaParabolaSlider
from m4.devices.reference_mirror import OpcUaReferenceMirror
from m4.ott_sim.fake_interferometer import FakeInterferometer
from m4.ott_sim.fake_accelerometers import FakeAccelerometers
from m4.ott_sim.fake_parabola_slider import FakeParabolaSlider
from m4.ott_sim.fake_reference_mirror import FakeReferenceMirror
from m4.devices.temperature_sensors import OpcUaTemperatureSensors
from m4.ott_sim.fake_temperature_sensors import FakeTemperatureSensors
from m4.devices.reference_mirror_slider import OpcUaReferenceMirrorSlider
from m4.ott_sim.fake_reference_mirror_slider import FakeReferenceMirrorSlider

def create_ott():
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
         - M4's Exapode
    interf: object
        The interferometer used for data acquisition.
    dm: object
        The deformable mirror, that is M4.
    """
    conf_obj = ufp.folders

    ###########
    # Sliders #
    ###########

    # Parabola Slider
    if conf_obj.simulated_parSlider is True:
        parabola_slider = FakeParabolaSlider()
    else:
        opcUa = OpcUaController()
        parabola_slider = OpcUaParabolaSlider(opcUa)

    # Reference Mirror Slider
    if conf_obj.simulated_rmSlider is True:
        reference_mirror_slider = FakeReferenceMirrorSlider()
    else:
        opcUa = OpcUaController()
        reference_mirror_slider = OpcUaReferenceMirrorSlider(opcUa)

    # Rotator
    if conf_obj.simulated_angleRotator is True:
        angle_rotator = FakeAngleRotator()
    else:
        opcUa = OpcUaController()
        angle_rotator = OpcUaAngleRotator(opcUa)

    #############
    # alignment #
    #############

    # Parabola
    if conf_obj.simulated_par is True:
        parab = FakeParabola()
    else:
        opcUa = OpcUaController()
        parab = OpcUaParabola(opcUa)

    # Reference Mirror
    if conf_obj.simulated_rm is True:
        reference_mirror = FakeReferenceMirror()
    else:
        opcUa = OpcUaController()
        reference_mirror = OpcUaReferenceMirror(opcUa)

    # DP
    if conf_obj.simulated_dp is True:
        dp = None
    else:
        dp = ZmqDpMotors()

    # M4
    if conf_obj.simulated_m4Exapode is True:
        m4 = FakeM4Exapode()
    else:
        opcUa = OpcUaController()
        m4 = OpcUaM4Exapode(opcUa)

    # DM
    if conf_obj.simulated_dm is True:
        dm = FakeM4DM()
    else:
        dm = None

    ##########    
    # Others #
    ##########

    # Temperature Sensors and Accelerometers
    if conf_obj.simulated_accelerometers is True:
        accelerometers = FakeAccelerometers()
    else:
        accelerometers = ZmqAccelerometers()

    if conf_obj.simulated_tempSensors is True:
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
    if conf_obj.simulated_interf is True:
        interf = FakeInterferometer(ott, dm)
    else:
        interf = I4d6110()

    if Sound.PLAY is True:
        playsound.playsound(os.path.join(Sound.AUDIO_FILE_PATH, "ott-ini.mp3"))
        playsound.playsound(os.path.join(Sound.AUDIO_FILE_PATH, "ott-conf.mp3"))

    return ott, interf, dm
