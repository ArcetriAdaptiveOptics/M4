import numpy as np
from m4.mini_OTT import timehistory as th
from m4.misc import par_filterCenter as pf
from m4.ground import zernike as zern
from importlib import reload
from m4.ground import geo
from matplotlib import pyplot as plt



######### 1st vs last set (to be specified) analysis
tn = '20230405_233640' #3x3
tn = '20230406_090110' #3x3 
tn = '20230406_140429' #3x3

zlist = [1,2,3,4]
aveimg, plist = pf.filterFrames(tn,rad1=0.1,rad2=0.02,thr=1.2)
dd = pf.diffAnalysis(tn,plist,bbreak=2,zlist=zlist)
cc, mat = zern.zernikeFit(dd,np.arange(11)+1)

plt.figure()
plt.imshow(dd, cmap='jet'); plt.colorbar();
plt.title(tn + ' RMS = %.3e m' %(dd.std()))

#in dd interessante provare a levare solo ptt





##########  OPDSeries analysis

tn = '20230405_233640' #3x3
tn = '20230406_090110' #3x3 
tn = '20230406_140429' #3x3




tn = '20230426_085332' # prima di sollevare la parabola
tn = '20230426_110517' # dopo aver sollevato la parabola   # infondo per cfr


thr = 1.2
zlist = [1,2,3,4,7,8,11]

par, plist = pf.filterFrames(tn,rad1=0.1,rad2=0.02,thr=thr)
parz = th.removeZernike(par,zlist)
plt.figure()
plt.clf(); plt.imshow(parz,cmap='jet'); plt.colorbar()
titolo = tn + ' RMS = %.3e m ' %(parz.std())
title(titolo)

#############################








######### 2x2 vs 3x3 kernel analysis #######

tn2x2 = '20230404_102354'
tn3x3 = '20230405_185019'

thr = 1.2
zlist = [1,2,3,4,7,8,11]

print('\n 2x2 kernel analysis')
par2x2, plist2x2 = pf.filterFrames(tn2x2,rad1=0.1,rad2=0.02,thr=thr)
par3x3, plist3x3 = pf.filterFrames(tn3x3,rad1=0.1,rad2=0.02,thr=thr)
print('\n 3x3 kernel analysis')
par3x3 = pf.filterFrames(tn3x3,rad1=0.1,rad2=0.02,thr=thr)


zlist = [1,2,3,4,7,8,11]
zlist = np.arange(11)+1

par2x2z = th.removeZernike(par2x2,zlist)
plt.figure(3)
plt.clf(); imshow(par2x2z,cmap='jet'); plt.colorbar()
titolo = tn2x2 + ' RMS = %.3e m ' %(par2x2z.std())
plt.title(titolo)

par3x3z = th.removeZernike(par3x3,zlist)
plt.figure(4)
plt.clf(); imshow(par3x3z,cmap='jet'); plt.colorbar()
titolo = tn3x3 + ' RMS = %.3e m ' %(par3x3z.std())
plt.title(titolo)


dd = par2x2z - par3x3z
rms_dd = dd.std()
print('\n \u0394RMS = %.3e m' %(rms_dd))

plt.figure(5)
plt.clf(); imshow(dd,cmap='jet',vmin=-2e-8,vmax=2e-8);plt.colorbar()
plt.title('\u0394 Average, RMS = %.3e m' %(rms_dd))


#####################################




#########     2023/04/26  SOLLEVAMENTO PARABOLA
######## analisi prima e dopo aver sollevato e riabbassato PAR
thr = 1.2
zlist = [1,2,3,4,7,8]
tn1 = '20230426_085332' # prima di sollevare la parabola
par, plist_1 = pf.filterFrames(tn1,rad1=0.1,rad2=0.02,thr=thr)
par_before = th.removeZernike(par,zlist)
tn2 = '20230426_110517' # dopo aver sollevato la parabola
par, plist_2 = pf.filterFrames(tn2,rad1=0.1,rad2=0.02,thr=thr)
par_after = th.removeZernike(par,zlist)

####### Creo maschera per il centro rad = .1 m
cir = geo.qpupil(np.invert(par_before.mask))
ps = pf.getPixelScale(tn1)
rad = .1*ps
inner_mask = geo.draw_mask(np.zeros(np.shape(par_before)),cir[0],cir[1],rad)

mmask = par_before.mask + inner_mask
par_before = np.ma.masked_array(par_before.data,mmask)
par_after = np.ma.masked_array(par_after.data,mmask)

dd = par_after - par_before



plt.figure(figsize=(18,5))
lim = 20e-9
plt.subplot(1,3,1); plt.subplot(1,3,1); plt.imshow(par_before,cmap='jet'); plt.colorbar();title('Par_before \n RMS = %.3e m' %(par_before.std()))
plt.subplot(1,3,2); plt.imshow(par_after,cmap='jet'); colorbar(); plt.title('Par_after  \n RMS = %.3e m' %(par_after.std()))
plt.subplot(1,3,3); plt.imshow(dd,vmin=-lim,vmax=lim,cmap='jet'); plt.colorbar(); title('\u0394Par RMS = %.3e m ' %(dd.std()))
plt.suptitle('PAR lifting tes')





