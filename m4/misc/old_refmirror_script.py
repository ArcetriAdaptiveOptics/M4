
import numpy as np
from matplotlib import pyplot as plt
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
fold='0000_0000_ref'

img=rmlib.patchdata(tn, fold)
cir=geo.qpupil(np.invert(img.mask))
refdiam = 0.1
pixscale = cir[2]*2/refdiam

ppos, fl = rmlib.fold2pos(tn)
i=0
img = rmlib.patchdata(tn, fl[i], thrsurf=2, thrpix=0.8)
iml = rmlib.getcube(rmlib.gimmethelist(tn, fl[i]))
imgvec, maskvec,ppos = rmlib.tndata(tn,thrsurf=2, thrpix=0.95)
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
plt.clf(); plt.imshow(pippo,vmin=-3e-7, vmax=3e-7); plt.colorbar()


fl = rmlib.gimmethelist(tn, fold)
cc = rmlib.getcube(fl)
cc1 = rmlib.mask_check(cc, 0.8)
cc2 = rmlib.std_check(cc1, 2)

#img = folddata(tn, fold, thrsurf=2, thrpix=0.8)



foldlist = os.listdir(base+tn)
pp,patchlist = rmlib.fold2pos(tn)
p0 = pp[0]
p1 = pp[1]

fullimgvec = rmlib.refmirr_fullframes(maskvec, imgvec,ppos, pixscale)


