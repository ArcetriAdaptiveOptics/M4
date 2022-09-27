'''
Authors
  - C. Selmi: written in 2020

HOW TO USE IT::
    from m4.configuration import start
    ott, interf = start.create_ott(config_file_path)
'''
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

def create_ott(config_file_name='/home/m4/git/M4/m4/configuration/towerConfig.yaml'):
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

    if conf_obj.simulated == 1:
        parabola_slider = FakeParabolaSlider()
        reference_mirror_slider = FakeReferenceMirrorSlider()
        angle_rotator = FakeAngleRotator()
        parab = FakeParabola()
        reference_mirror = FakeReferenceMirror()
        m4 = FakeM4Exapode()
        temperature_sensor = FakeTemperatureSensors()
        accelerometers = FakeAccelerometers()
        interf = FakeInterferometer()
        dm = FakeM4DM()
    else:
        opcUa = OpcUaController()
        parabola_slider = OpcUaParabolaSlider(opcUa)
        reference_mirror_slider = OpcUaReferenceMirrorSlider(opcUa)
        angle_rotator = OpcUaAngleRotator(opcUa)
        parab = OpcUaParabola(opcUa)
        reference_mirror = OpcUaReferenceMirror(opcUa)
        m4 = OpcUaM4Exapode(opcUa)
        temperature_sensor = OpcUaTemperatureSensors(opcUa)
        accelerometers = ZmqAccelerometers()
        interf = I4d6110()
        dm = None


    ott = OTT(parabola_slider, reference_mirror_slider, angle_rotator,
              parab, reference_mirror, m4, temperature_sensor, accelerometers)
    if conf_obj.simulated == 1:
        interf.set_ott(ott)
        interf.set_dm(dm)

    return ott, interf, dm
