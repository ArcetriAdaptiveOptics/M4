conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
from m4.configuration import start
import numpy as np
import matplotlib.pyplot as plt
ott, interf, dm = start.create_ott(conf)
from m4.configuration.ott_parameters import Interferometer
from m4.mini_OTT import timehistory as th

print("Frequency acquired: " , Interferometer.BURST_FREQ)

from m4 import noise

comp_stf=False
comp_ast=True

tnlist=['20230126_125637','20230126_130004','20230126_131049','20230126_132012', '20230126_160816']
tn = tnlist[4]
if comp_stf:
    #structure function
    for tn in tnlist:
        dfpath =  '/home/m4/4d/M4/Produced/'+tn
        print(dfpath)
        tau_vector = np.arange(1,100,2)
        noise.convection_noise(dfpath, tau_vector)


if comp_ast:
    zv = []
    q = []
    #astigmatism analysis
    for tn in tnlist:
        fl = th.fileList(tn)
        img = th.averageFrames(0,499,fl)
        zz, mat =th.zernike.zernikeFit(img, [1,2,3,4,5,6])
        img = th.removeZernike(img)
        q.append(img)
        print(tn, zz)
        zv.append(zz)

    zv = np.array(zv)
        

