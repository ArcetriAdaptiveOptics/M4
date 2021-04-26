'''
Authors
  - C. Selmi: written in 2020
'''
from m4.configuration.create_ott import OTT
from m4.ground.interface_4D import comm4d
from m4.configuration import config
from m4.ground.opc_ua_controller import OpcUaController
from m4.ott_sim.fake_parabola_slider import FakeParabolaSlider
from m4.ott_sim.fake_reference_mirror_slider import FakeReferenceMirrorSlider
from m4.ott_sim.fake_angle_rotator import FakeAngleRotator
from m4.devices.parabola_slider import OpcUaParabolaSlider
from m4.devices.reference_mirror_slider import OpcUaReferenceMirrorSlider
from m4.devices.angle_rotator import OpcUaAngleRotator
from m4.devices.parabola import OpcUaParabola

def create_ott():
    ''' Function for the ott creation

    Returns
    -------
    ott: object
        tower
    interf: object
        interferometer
    '''
    if config.simulated == 1:
        parabola_slider = FakeParabolaSlider()
        reference_mirror_slider = FakeReferenceMirrorSlider()
        angle_rotator = FakeAngleRotator()
        parabola = FakeParabolaSlider()
    else:
        opcUa = OpcUaController()
        parabola_slider = OpcUaParabolaSlider(opcUa)
        reference_mirror_slider = OpcUaReferenceMirrorSlider(opcUa)
        angle_rotator = OpcUaAngleRotator(opcUa)
        parabola = OpcUaParabola(opcUa)

    ott = OTT(parabola_slider, reference_mirror_slider, angle_rotator, parabola)
    interf = comm4d()
    return ott, interf
