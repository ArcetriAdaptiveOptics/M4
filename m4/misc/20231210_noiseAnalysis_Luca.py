'''
This code is meant to perform the noise analysis on the burst acuisitions.

The 3 different plots that produce are:

the absolute noise error
DESCRIVO

the vibration error
DESCRIVO

the piston noise convection error
DESCRIVO

'''
import os
import numpy as np
from matplotlib.pyplot import *
from m4.configuration import start
from m4.mini_OTT import timehistory as th
from m4.ground import zernike
from m4.ground import timestamp
import time as tt
from m4 import noise

##

tn=[]
tn.append('20231011_144711')
tn.append('20231007_113703')
tn.append('20231006_232601')
tn.append('20231005_212609')
tn.append('20231028_084255')

#

tau_vector = np.arange(1,100,2)
path_series = '/mnt/m4storage/Data/M4Data/OPTData/OPDImages/'

close('all')

for jj in tn:
    dfpath=os.path.join(path_series,jj)
    noise.convection_noise(dfpath, tau_vector)
    

template=np.array([3,11,25,37,51])

path_series = '/mnt/m4storage/Data/M4Data/OPTData/OPDImages/'

close('all')

for jj in tn:
    dfpath=os.path.join(path_series,jj)
    noise.noise_vibrations(dfpath,template,0)
    
close('all')
