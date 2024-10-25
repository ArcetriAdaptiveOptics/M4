'''
This code is meant to perform the noise analysis on the burst acuisitions.


'''
# import os
# import numpy as np
# from matplotlib.pyplot import *
# from m4.configuration import start
# from m4.mini_OTT import timehistory as th
# from m4.ground import zernike
# from m4.ground import timestamp
# import time as tt
from m4 import noise

##

tn=[]
tn.append(['20240521_152923', 55.5]) # damper off
tn.append(['20240521_153943', 440.]) # damper off
tn.append(['20240521_160030', 55.5]) # damper on
tn.append(['20240521_155828', 440.]) # damper on

# tn.append('20231007_113703')
# tn.append('20231006_232601')
# tn.append('20231005_212609')
# tn.append('20231028_084255')

#
## rumore di convezione
tau_vector = np.arange(1,100,2)
path_series = '/mnt/m4storage/Data/M4Data/OPTData/OPDImages/'

plt.close('all')

for jj in tn:
    dfpath=os.path.join(path_series,jj[0])
    noise.convection_noise(dfpath, tau_vector, freq=jj[1])

## rumori di vibrazione ????
template=np.array([3,11,25,37,51])

path_series = '/mnt/m4storage/Data/M4Data/OPTData/OPDImages/'

plt.close('all')

for jj in tn:
    dfpath = th.foldname.OPD_IMAGES_ROOT_FOLDER+'/'+jj[0]+'/'
    noise.noise_vibrations(dfpath,template,0)

##
