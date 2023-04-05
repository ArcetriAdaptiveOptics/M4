
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
'''
analysis template
the tracknum contains an entire, complete and consistent dataset.
the dataset is organized as:
tn/000-000-ref0  tn/000-000-ref1
tn/000-000  tn/000-001  tn/000-002
tn/001-000  tn/001-001  tn/001-002

the following inputs shall be passed: pixelscale; coordinates of the center of frames

algorithm:
tn is passed
folders are identified
ref images are prepared - ref mask are identified
patch masks are extracted and cropped to ref mask
global mask is created from patch masks
patch images are extracted, filtered and averaged
global image is created

'''

base = '/mnt/m4storage/Data/M4Data/OPTData/RefMirror/'
tn = '20230315_155719'

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

ave1 = np.ma.mean(imgvec[0:2],0)
ave2 = np.ma.mean(imgvec[3:5],0)
ave = (ave1+ave2)/2


for ii in range(len(imgvec)):
    imgvec[ii] = imgvec[ii] - ave
#imgvec = imgvec2



#### subiap img check
seq = [2,3,0,1,4,5]
plt.figure(figsize=(20,20))
vlim = .3e-7
for ii in range(len(seq)):
    plt.subplot(2,3,ii+1)
    plt.imshow(imgvec[seq[ii]],vmin=-vlim, vmax=vlim,cmap='jet')
    plt.colorbar()
    plt.title(ii)



##### STITCH
ppos=np.array([[2,0],[3,0],[0,0],[1,0],[4,0],[5,0]])
ppos = ppos*0.1
fullmask = rmlib.refmirr_fullmask(maskvec, ppos, pixscale)
fullimgvec = rmlib.refmirr_fullframes(maskvec, imgvec,ppos, pixscale)

seq = [2,3,0,1,4,5]

iv = []
for i in seq:
    iv.append(fullimgvec[i])

fullframe = rmlib.stitch_process2(iv, fullmask,zlist2fit=[1,2,3],zlist2adjust=[1,2,3])


pippo = th.removeZernike(fullframe,[1,2,3])
fig=plt.figure(figsize=(5,10)); plt.imshow(pippo,cmap='jet',vmin=-10e-8, vmax=10e-8); plt.colorbar(); plt.title('RMS = %.3e nm' %(np.std(pippo)))








