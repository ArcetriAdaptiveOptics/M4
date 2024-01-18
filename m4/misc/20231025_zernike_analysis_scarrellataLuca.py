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


##  define the tracking numbers 
tn0='20231010_130115'
tn1 ='20231010_190535'

tnlist = th.tnscan(tn0, tn1)
ntn = len(tnlist)

## open average if already produced
qc = []
for i in tnlist:
    # q = th.openAverage(i)
    # #plt.figure();plt.imshow(q);plt.clim([-10e-8, 10e-8]);plt.colorbar();plt.show()
    # q = th.removeZernike(q,[1,2,3,4])
    q=imgreg.load_ott(i, zlist=[1,2,3])
    qc.append(q) 

## Here we compute the zernike coefficient for focus and astigmatism and coma for each sub apertures
#plt.close('all')

tnpar  = '20231016_124531' #old '20231011_150952'
px_ott = 0.00076
f0 = 0
f1 = 40
par_remapped = imgreg.load_registeredPar(tnpar)
par_filtered = th.comp_filtered_image(par_remapped,  d=px_ott, verbose=True, disp=False, freq2filter=(f0,f1))

par= par_remapped

mm = geo.draw_mask(np.ones(np.shape(par.data)), np.shape(par)[0]/2, np.shape(par)[1]/2, 800, out=1)
par.mask[mm!=0]=1 
# plt.figure();plt.imshow(par); plt.show()
# 
# i=0
# mask=qc[i].mask
# par.mask=mask
# plt.figure(); plt.imshow(par); plt.show()
# cc, mm = zernike.zernikeFit(par, [1, 2, 3, 4, 5, 6, 7, 8])
# print(cc[6:8])
# np.std(par)

coeff_ott=[]
coeff_par=[]
coeff_res=[]
k1=0;k2=1
for i in  np.arange(k1,ntn-k2,1):
    
    #scelgo se usare la parabola rimappata o filtrata
    par = imgreg.image_remask(par, qc[i])
    res = qc[i]-2*par
      
    cc,mat=zernike.zernikeFit(qc[i],[1,2,3,4,5,6,7,8])
    coeff_ott.append(cc)

    cc,mat=zernike.zernikeFit(2*par,[1,2,3,4,5,6,7,8])
    coeff_par.append(cc)
    
    cc,mat=zernike.zernikeFit(res,[1,2,3,4,5,6,7,8])
    coeff_res.append(cc)
    
    # faccio il fit degli zernike sulla pupilla grossa invece che su quella piccola. NON FUNZIONA?
    # cc,mat=zernike.zernikeFitAuxmask(res, par_remapped.mask, [1,2,3,4,5,6,7,8])
    # coeff_aux.append(cc)
    

coeff_ott=np.array(coeff_ott)
coeff_par=np.array(coeff_par)
coeff_res=np.array(coeff_res)
fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(16,5))
ax1.plot(coeff_ott[:,3]);ax1.plot(coeff_ott[:,4]);ax1.plot(coeff_ott[:,5]);ax1.plot(coeff_ott[:,6]);ax1.plot(coeff_ott[:,7])#; ax1.set_ylim([-5e-8,5e-8])
ax2.plot(coeff_par[:,3]);ax2.plot(coeff_par[:,4]);ax2.plot(coeff_par[:,5]);ax2.plot(coeff_par[:,6]);ax2.plot(coeff_par[:,7])#; ax2.set_ylim([-5e-8,5e-8]) 
ax3.plot(coeff_res[:,3]);ax3.plot(coeff_res[:,4]);ax3.plot(coeff_res[:,5]);ax3.plot(coeff_res[:,6]);ax3.plot(coeff_res[:,7])#; ax3.set_ylim([-5e-8,5e-8])
ax1.set_title('ott'); ax2.set_title('2*par'); ax3.set_title('res=ott-2*par')
ax1.legend(['z4','z5','z6','z7','z8']);ax2.legend(['z4','z5','z6','z7','z8']);ax3.legend(['z4','z5','z6','z7','z8']); fig.suptitle('Zernike coefficient for different positions in the tower')
plt.show()
