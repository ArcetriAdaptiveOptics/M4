'''
Authors
  - C. Selmi:  written in 2021
'''

def ciccio():
    from m4.configuration.create_ott import OTT
    ott = OTT()

    angle = ott.angle()
    rslide = ott.rslide()
    slide = ott.slide()
    par = ott.parab()
    rm = ott.refflat()
    temperature = ott.temperature()
