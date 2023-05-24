
import numpy as np
from importlib import reload
#import sys
from astropy.io import fits as pyfits
#from subapcalib import sasimulator as sasim
#from subapcalib import safunct as sa
from m4.ground import zernike as zern
from m4.ground import geo
from m4.mini_OTT import timehistory as th
from m4.misc import refmirror_lib as rmlib
import os
from matplotlib import pyplot as plt

base = '/mnt/m4storage/Data/M4Data/OPTData/RefMirror/'
#tn = '20230328_163013'
tn = '20230331_124514'
step = .125 #0.06  #[m]

######  GET REF DATA ########  Check rmlib.py for cavity definition
pixscale = rmlib.getrefdata(tn)

###### GET DATA FROM MEASUREMENT #####
ppos, fl = rmlib.fold2pos(tn)
imgvec0, maskvec0,ppos = rmlib.tndata(tn,thrsurf=2, thrpix=0.95,rebinfactor=2)

###### ERODE MASK AND REMOVE PTT  (problemi nell'erosione, maskera ovale)
imgvec, maskvec = rmlib.erodemask(maskvec0,imgvec0,dpix=2)

imgvec2 = []
for ii in range(len(imgvec)):
     aa = th.removeZernike(imgvec[ii],[1,2,3])
     imgvec2.append(aa)
imgvec = imgvec2

##### AVERAGE ########
ave = np.ma.mean(imgvec,axis=0)



for ii in range(len(imgvec)):
    imgvec[ii] = imgvec[ii] - ave
#imgvec = imgvec2



#### subiap img check
plt.figure(figsize=(20,20))
vlim = 1e-7
for ii in range(len(imgvec)):
    plt.subplot(3,7,ii+1)
    plt.imshow(imgvec[ii],cmap='jet')#vmin=-vlim, vmax=vlim)
    plt.colorbar()
    plt.title(ii)



##### STITCH
#ppos=np.array([[2,0],[3,0],[0,0],[1,0],[4,0],[5,0]])
#ippos = ppos*step
ppos = np.flip(ppos,1)
aa = np.array([4, 4]) - ppos
ppos = aa*step


fullmask = rmlib.refmirr_fullmask(maskvec, ppos, pixscale)
fullimgvec = rmlib.refmirr_fullframes(maskvec, imgvec,ppos, pixscale)

iv = fullimgvec

zlist = np.arange(4)+1

fullframe = rmlib.stitch_process2(iv, fullmask,zlist2fit=zlist,zlist2adjust=zlist)

pippo = th.removeZernike(fullframe,zlist)
fig=plt.figure(figsize=(10,10)); plt.imshow(pippo,cmap='jet'); plt.colorbar(); plt.title('RMS = %.3e m' %(np.std(pippo)))


#------ analysis of in-out with ref and center
tn = '20230418_121136'
imgv, maskv, refv, maskrefv =  rmlib.tndataall(tn, thrsurf=2, thrpix=0.8, rebinfactor=8)
zr = []
zm = []
for i in imgv:
    cc, m= zern.zernikeFit(i, [1,2,3,4,5,6,7,8,9,10,11])
    zm.append(cc)
for i in refv:
    cc, m= zern.zernikeFit(i, [1,2,3,4,5,6,7,8,9,10,11])
    zr.append(cc)
zr = np.array(zr)
zm = np.array(zm)
plot(zm[:,3], 'o')
plot(zr[:,3],'x')
refflat_power=11.3e-9
refmirr = np.array([ -2.11e-09,   1.123e-06,   4.110e-07,   1.12989e-08,   5.45e-10,  -1.25282e-09,  -1.15052e-10,  -1.64e-09 , -1.120e-10,   3.484e-10 , -5.058e-09])
cc = []
for i in range(len(refmirr)):
    cc.append( np.mean(zm[:,i])-np.mean(zr[:,i])+refmirr[i])

print(cc)



