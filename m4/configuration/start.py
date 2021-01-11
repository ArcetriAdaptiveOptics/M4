'''
Authors
  - C. Selmi: written in 2020
'''
from m4.configuration.create_ott import OTT
from m4.ground.interface_4D import comm4d

def create_ott():
    ''' Function for the ott creation
    '''
    ott = OTT()
    interf = comm4d()
    return ott, interf
