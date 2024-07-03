'''
Authors
  - C. Selmi: written in 2020
              modified in 2022

HOW TO USE IT::
    from m4.configuration import start
    ott, interf, dm = start.create_ott(config_file_path)
'''
import os
from m4.configuration.create_ott import OTT
from m4.configuration.config_reader import configuration_path
from m4.devices.opc_ua_controller import OpcUaController
from m4.ott_sim.fake_parabola_slider import FakeParabolaSlider
from m4.ott_sim.fake_reference_mirror_slider import FakeReferenceMirrorSlider
from m4.ott_sim.fake_angle_rotator import FakeAngleRotator
from m4.ott_sim.fake_parabola import FakeParabola
from m4.ott_sim.fake_reference_mirror import FakeReferenceMirror
from m4.ott_sim.fake_m4_exapode import FakeM4Exapode
from m4.ott_sim.fake_temperature_sensors import FakeTemperatureSensors
from m4.ott_sim.fake_interferometer import FakeInterferometer
from m4.ott_sim.fake_accelerometers import FakeAccelerometers
from m4.ott_sim.fake_deformable_mirror import FakeM4DM
from m4.devices.parabola_slider import OpcUaParabolaSlider
from m4.devices.reference_mirror_slider import OpcUaReferenceMirrorSlider
from m4.devices.angle_rotator import OpcUaAngleRotator
from m4.devices.parabola import OpcUaParabola
from m4.devices.reference_mirror import OpcUaReferenceMirror
from m4.devices.m4_exapode import OpcUaM4Exapode
from m4.devices.temperature_sensors import OpcUaTemperatureSensors
from m4.devices.accelerometers import ZmqAccelerometers
from m4.devices.interferometer import I4d6110

from m4.configuration.config_uploader import config_rewriter
from m4.configuration.ott_parameters import Sound
import playsound

def create_ott(config_file_name='/home/m4/git/M4/m4/configuration/myConfig.yaml'):
    ''' Function for the ott creation

    Parameters
    ---------
    config_file_name: string
        configuration path to use

    Returns
    -------
    ott: object
        tower
    interf: object
        interferometer
    '''
    conf_obj = configuration_path(config_file_name)
    cr = config_rewriter(conf_obj)
    cr.upload()

    if conf_obj.simulated_accelerometers is True:
        accelerometers = FakeAccelerometers()
    else:
        accelerometers = ZmqAccelerometers()
    if conf_obj.simulated_angleRotator is True:
        angle_rotator = FakeAngleRotator()
    else:
        opcUa = OpcUaController()
        angle_rotator = OpcUaAngleRotator(opcUa)
    if conf_obj.simulated_m4Exapode is True:
        m4 = FakeM4Exapode()
    else:
        opcUa = OpcUaController()
        m4 = OpcUaM4Exapode(opcUa)
    if conf_obj.simulated_parSlider is True:
        parabola_slider = FakeParabolaSlider()
    else:
        opcUa = OpcUaController()
        parabola_slider = OpcUaParabolaSlider(opcUa)
    if conf_obj.simulated_par is True:
        parab = FakeParabola()
    else:
        opcUa = OpcUaController()
        parab = OpcUaParabola(opcUa)
    if conf_obj.simulated_rmSlider is True:
        reference_mirror_slider = FakeReferenceMirrorSlider()
    else:
        opcUa = OpcUaController()
        reference_mirror_slider = OpcUaReferenceMirrorSlider(opcUa)
    if conf_obj.simulated_rm is True:
        reference_mirror = FakeReferenceMirror()
    else:
        opcUa = OpcUaController()
        reference_mirror = OpcUaReferenceMirror(opcUa)
    if conf_obj.simulated_tempSensors is True: 
        temperature_sensor = FakeTemperatureSensors()
    else:
        opcUa = OpcUaController()
        temperature_sensor = OpcUaTemperatureSensors(opcUa)


    if conf_obj.simulated_interf is True:
        interf = FakeInterferometer()
    else:
        interf = I4d6110()

    if conf_obj.simulated_dm is True:
        dm = FakeM4DM()
    else:
        dm = None

    ott = OTT(parabola_slider, reference_mirror_slider, angle_rotator,
              parab, reference_mirror, m4, temperature_sensor, accelerometers)

    if Sound.PLAY is True:
        playsound.playsound(os.path.join(Sound.AUDIO_FILE_PATH,'ott-ini.mp3'))

    if conf_obj.simulated_interf is True:
        interf.set_ott(ott)
        interf.set_dm(dm)
    if Sound.PLAY is True:
        playsound.playsound(os.path.join(Sound.AUDIO_FILE_PATH, 'ott-conf.mp3'))
    return ott, interf, dm
