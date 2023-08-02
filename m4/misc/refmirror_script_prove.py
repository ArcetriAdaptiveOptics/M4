
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
aa = np.where(ppos[:,0]==1)
ppos = ppos[aa]
foldref = '0001_0000_ref0'
refimg = rmlib.patchdata(tn,foldref, thrsurf=2, thrpix=0.8, rebinfactor=8)
refimg = th.removeZernike(refimg,[1,2,3])
imshow(refimg,cmap='jet')

imgvec0, maskvec0,ppos = rmlib.tndata(tn,thrsurf=2, thrpix=0.95)
ppos = ppos[aa]
img = []
mask = []
for ii in aa[0]:
     img.append(imgvec0[int(ii)])
     mask.append(maskvec0[int(ii)])


###### ERODE MASK AND REMOVE PTT
imgvec, maskvec = rmlib.erodemask(mask,img,dpix=2)

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
figure(figsize=(20,20))
vlim = 1e-7
for ii in range(len(imgvec)):
    subplot(3,7,ii+1)
    imshow(imgvec[ii],cmap='jet')#vmin=-vlim, vmax=vlim)
    colorbar()
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

fullframe = rmlib.stitch_process2(iv, fullmask,zlist2fit=[1,2,3],zlist2adjust=[1,2,3])


pippo = th.removeZernike(fullframe,[1,2,3])
fig=plt.figure(figsize=(10,10)); plt.imshow(pippo,cmap='jet',vmin=-10e-8, vmax=10e-8); plt.colorbar(); plt.title('RMS = %.3e nm' %(np.std(pippo)))








