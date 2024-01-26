import numpy as np
from m4.mini_OTT import timehistory as th
from m4.misc import par_filterCenter as pf
from m4.ground import zernike as zern
from importlib import reload
from m4 import requirements_checker as rc
from m4.analyzers import requirement_analyzer as ra
from m4.ground import geo
from m4.misc import refmirror_lib as rmlib




#### take 1 image to analyze ####
tn = '20230405_233640' #3x3
aveimg, plist = pf.filterFrames(tn,rad1=0.1,rad2=0.02,thr=1.2)
#### pixel scale
pscale = pf.getPixelScale(tn)
############
##### data for ReqCheck
step = 2000000
n_patches = None
zernlist = [1,2,3]




######## Simple req check PTT removed ######
image = th.removeZernike(aveimg,zernlist)
slope = ra.test242(image, pscale)
diff_piston = ra.diffPiston(image)
roc = ra.test283(image, pscale, step)
rms31 = ra.test243(image, 0.015, pscale, step, n_patches)
rms500 = ra.test243(image, 0.1, pscale, step, n_patches)
############################################
###### Req Check w/ rebinning PTT removed ##
factor = 2
image = rmlib.rebin(aveimg,factor)
image = th.removeZernike(image,zernlist)
slope = ra.test242(image, pscale/factor)
diff_piston = ra.diffPiston(image)
roc = ra.test283(image, pscale/factor, step)
rms31 = ra.test243(image, 0.015, pscale/factor, step, n_patches)
rms500 = ra.test243(image, 0.1, pscale/factor, step, n_patches)
###########################################
##### Req check after nxn shift and subtraction
p_shift = 1
image = th.removeZernike(aveimg,zernlist)
shift = np.roll(image,(p_shift,p_shift),(0,1))
diff_ima = image - shift
slope = ra.test242(diff_ima, pscale)
diff_piston = ra.diffPiston(diff_ima)
roc = ra.test283(diff_ima, pscale, step)
rms31 = ra.test243(diff_ima, 0.015, pscale, step, n_patches)
rms500 = ra.test243(diff_ima, 0.1, pscale, step, n_patches)
##########################################
#### Req check shift + rebin #############
p_shift = 1
factor = 2
image = th.removeZernike(aveimg,zernlist)
shift = np.roll(image,(p_shift,p_shift),(0,1))
diff_ima = image - shift
reb_ima = rmlib.rebin(diff_ima,factor)
slope = ra.test242(reb_ima, pscale/factor)
diff_piston = ra.diffPiston(reb_ima)
roc = ra.test283(reb_ima, pscale/factor, step)
rms31 = ra.test243(reb_ima, 0.015, pscale/factor, step, n_patches)
rms500 = ra.test243(reb_ima, 0.1, pscale/factor, step, n_patches)
########################################











