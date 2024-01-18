
from m4.configuration import start
from m4.mini_OTT import timehistory as th
from m4.ground import zernike
from m4.ground import timestamp
import time as tt
from m4 import noise
import os
import numpy as np
from matplotlib import pyplot as plt
from m4.ground import geo
from importlib import reload
from m4.misc import image_registration_lib as imgreg
from astropy.io import fits as pyfits


##
tn='20231019_002121'


plt.close('all')


fl = th.fileList(tn)
nf = len(fl)
q0=th.averageFrames(0,999,fl)

plt.figure(); plt.imshow(th.removeZernike(q0,[1,2,3,4,7,8])); plt.show()

##
rr = []
cc = []
kk=50
for i in np.arange(0,int(999/kk),1):
    print(i)
    #qi = th.removeZernike(th.frame(i,fl)-q0,[1,2,3,4,7,8])
    ave=th.averageFrames(i*kk,i*kk+kk,fl)
    c,mat=zernike.zernikeFit(ave,[1,2,3,4,5,6])
    cc.append(c)
    
    #rr.append(qi.std())
    
    ##

sec=30; min=sec/60; hour=min/60
xx=np.arange(0,len(cc),1)*hour*kk
cc=np.array(cc)

plt.figure()
plt.plot(xx,cc[:,3]);plt.plot(xx,cc[:,4]);plt.plot(xx,cc[:,5])
plt.legend(['z4','z5','z6'])
plt.title("zernike stability, tn="+tn)
plt.show()    

plt.figure()
plt.plot(xx,cc[:,3]-np.mean(cc[:,3]));plt.plot(xx,cc[:,4]-np.mean(cc[:,4]));plt.plot(xx,cc[:,5]-np.mean(cc[:,5]))
plt.xlabel('hours'); plt.ylabel('m')
plt.legend(['z4-meanz4','z5-meanz5','z6-meanz6'])
plt.title("zernike stability, tn="+tn)
plt.show()    