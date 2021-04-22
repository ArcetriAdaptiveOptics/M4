'''
Authors
  - C. Selmi: written in 2020
'''
from m4.configuration.create_ott import OTT
from m4.ground.interface_4D import comm4d
from m4.ott_sim.fake_parabola_slider import FakeParabolaSlider
from m4.ground.opc_ua_controller import OpcUaController
from m4.devices.parabola_slider import OpcUaParabolaSlider


def create_ott():
    ''' Function for the ott creation

    Returns
    -------
    ott: object
        tower
    interf: object
        interferometer
    '''
    if conf.simulated == 1:
        parabola_slider = FakeParabolaSlider()
    else:
        opcUa = OpcUaController()
        parabola_slider = OpcUaParabolaSlider(opcUa)

    ott = OTT(parabola_slider)
    interf = comm4d()
    return ott, interf
