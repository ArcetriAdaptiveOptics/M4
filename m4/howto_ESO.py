#from m4ott: calpyott
import m4
from m4.devices import opt_beam

ott,_,interf=m4.create_ott()

par      = opt_beam.Parabola(ott)
refmirr = opt_beam.ReferenceMirror(ott)






#from micws
##    calpy -f /data/Arcetri/Data/M4Data/SysConfig/configuration.yam
import numpy as np
from opticalib import PhaseCam, AdOpticaDm
from opticalib.dmutils import iff_module as ifm, iff_processing as ifp, iff_acquisition_preparation as ifa
from opticalib.dmutils.flattening import Flattening
from opticalib import load_fits, save_fits
from matplotlib.pyplot import *
ion()
