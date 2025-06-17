### restore the following lines to implement M4 init
'''
from opticalib import PhaseCam, AdOpticaDM
from m4.configuration import ott

interf = PhaseCam('6110')
dm = AdOpticaDM()

ott = ott.create_ott()
'''
