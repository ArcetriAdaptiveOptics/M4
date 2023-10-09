
import numpy as np
from importlib import reload
from m4.ground import zernike as zern
from m4.ground import geo
from m4.mini_OTT import timehistory as th
import os
from matplotlib import pyplot as plt
from m4.misc import par_filterCenter as pf




## Tn parabola w/cgh
tn = '20230501_152955'
fl=th.fileList(tn)
img1=th.averageFrames(0,300,fl)
img1=th.removeZernike(img1,[1,2,3,4])
img1 = np.flip(img1,axis=1)
ps = pf.getPixelScale(tn)   #pixelscale pix/m

## new mask definition
cir = geo.qpupil(np.invert(img1.mask))
rad = ps*.30
cc = ((cir[0])/2-17 ,cir[1])
newmask = geo.draw_mask(img1.data,cc[0],cc[1],rad)
newmask = -newmask+1
img2 = np.ma.masked_array(img1.data,newmask)
img2 = th.removeZernike(img2,[1,2,3,4])
plt.imshow(img1.mask)
plt.imshow(img2)
cc_par, mat_par = zern.zernikeFit(img2, np.arange(36)+1)



#############################
##### CHECK REQUIREMENTS ####
#############################

from m4.analyzers import requirement_analyzer as ra

## data for req check
step = 10000
n_patches = None
zernlist = [1,2,3,4]


#########
## rm al bordo di par

tn = '20230704_194022'
fl = th.fileList(tn)
img1 = th.averageFrames(0,100,fl)
img1 = th.removeZernike(img1, [1,2,3,4])
cir = geo.qpupil(np.invert(img1.mask))
### sottraggo i modi di par
cc, mat = zern.zernikeFit(img1, np.arange(36)+1)
surf = zern.zernikeSurface(img1, cc_par, mat)
dd = img1 - surf
dd = img1

### serve fare crop dell'immagine)
pp = cir[2]+20
img2 = img1[int(cir[0]-pp):int(cir[0]+pp),int(cir[1]-pp):int(cir[1]+pp)]
img2 = th.removeZernike(img2,zernlist)


## requirements script
ps = cir[2]*2/0.6  #pscale pix/m
slope = ra.test242(img2, ps)
diff_piston = ra.diffPiston(img2)
roc = ra.test283(img2,ps,step)
rms31 = ra.test243(img2, 0.015, ps, step, n_patches)
rms500 = ra.test243(img2, 0.1, ps, step, n_patches)
print('############### RESULTS ########### \nslope = %.3e asec \ndiff Piston = %.3e nm \nRoC = %.3e m \nrms31 = %.3e nm \nrms500 = %.3e nm' %(slope, diff_piston, roc, rms31, rms500))


###########################
## centro

tn = '20230706_192127'

fl = th.fileList(tn)
img1 = th.averageFrames(0,100,fl)
img1 = th.removeZernike(img1, [1,2,3,4])
cir = geo.qpupil(np.invert(img1.mask))
### serve fare crop dell'immagine)
pp = cir[2]+20
img2 = img1[int(cir[0]-pp):int(cir[0]+pp),int(cir[1]-pp):int(cir[1]+pp)]
img2 = th.removeZernike(img2,zernlist)
## requirements script
ps = cir[2]*2/0.5  #pscale pix/m
slope = ra.test242(img2, ps)
diff_piston = ra.diffPiston(img2)
roc = ra.test283(img2,ps,step)
rms31 = ra.test243(img2, 0.015, ps, step, n_patches)
rms500 = ra.test243(img2, 0.1, ps, step, n_patches)
print('############### RESULTS ########### \nslope = %.3e asec \ndiff Piston = %.3e nm \nRoC = %.3e m \nrms31 = %.3e nm \nrms500 = %.3e nm' %(slope, diff_piston, roc, rms31, rms500))


##########################
## fine corsa

tn = '20230707_104616'

fl = th.fileList(tn)
img1 = th.averageFrames(0,9,fl)
img1 = th.removeZernike(img1, [1,2,3,4])
cir = geo.qpupil(np.invert(img1.mask))
### serve fare crop dell'immagine)
pp = cir[2]+20
img2 = img1[int(cir[0]-pp):int(cir[0]+pp),int(cir[1]-pp):int(cir[1]+pp)]















