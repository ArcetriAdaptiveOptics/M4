
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

base = '/mnt/m4storage/Data/M4Data/OPTData/RefMirror/'
tn = '20230329_091725'
tn = '20230328_163013'
step = 0.06  #[m]

######  GET REF DATA ########  Check rmlib.py for cavity definition
pixscale = rmlib.getrefdata(tn)


###### GET DATA FROM MEASUREMENT #####
ppos, fl = rmlib.fold2pos(tn)
imgvec0, maskvec0,ppos = rmlib.tndata(tn,thrsurf=2, thrpix=0.95)

###### ERODE MASK AND REMOVE PTT
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
#figure(figsize=(20,20))
#vlim = .3e-7
#for ii in range(len(seq)):
#    subplot(2,3,ii+1)
#    imshow(imgvec[seq[ii]],vmin=-vlim, vmax=vlim,cmap='jet')
#    colorbar()
#    plt.title(ii)



##### STITCH
#ppos=np.array([[2,0],[3,0],[0,0],[1,0],[4,0],[5,0]])
#ippos = ppos*step
ppos = flip(ppos,1)
aa = np.array([10, 8]) - ppos
ppos = aa*step


fullmask = rmlib.refmirr_fullmask(maskvec, ppos, pixscale)
fullimgvec = rmlib.refmirr_fullframes(maskvec, imgvec,ppos, pixscale)

iv = fullimgvec

fullframe = rmlib.stitch_process2(iv, fullmask,zlist2fit=[1,2,3],zlist2adjust=[1,2,3])


pippo = th.removeZernike(fullframe,[1,2,3])
fig=figure(figsize=(5,10)); imshow(pippo,cmap='jet',vmin=-10e-8, vmax=10e-8); colorbar(); title('RMS = %.3e nm' %(std(pippo)))








