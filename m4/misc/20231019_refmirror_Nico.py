
from m4.configuration import start
from m4.mini_OTT import timehistory as th
from m4.ground import zernike
from m4.ground import timestamp
import time as tt
from m4 import noise
import os
import numpy as np
from matplotlib.pyplot import *
from m4.ground import geo
from importlib import reload
from m4.misc import image_registration_lib as imgreg
from astropy.io import fits as pyfits

##  define the tracking numbers 
tn0='20231010_130115'       # muovendo RM
tn1 ='20231010_190535'

tnlist = th.tnscan(tn0, tn1)
ntn = len(tnlist)

qc = []
for ii in tnlist:
    # q = th.openAverage(i)
    # #plt.figure();plt.imshow(q);plt.clim([-10e-8, 10e-8]);plt.colorbar();plt.show()
    # q = th.removeZernike(q,[1,2,3,4])
    q=imgreg.load_ott(ii, zlist=[1,2,3,4])
    cir = geo.qpupil(-1*q.mask+1)
    mm = geo.draw_mask(np.ones(np.shape(q.data)),cir[0],cir[1],cir[2]-3,out=1)
    q = np.ma.masked_array(q.data,mm)
    qc.append(q)
qc = np.ma.masked_array(qc)

ave = np.ma.average(qc,axis=0)
sum_mask = np.sum(-1*qc.mask+1,axis=0)
id_max = np.where(sum_mask==np.max(sum_mask))
mask = np.ones(np.shape(ave.data))
mask[id_max] = 0
eye_ave = np.ma.masked_array(ave.data, mask)

# Load par data, remap and filter
tnpar  = '20231016_124531' #old '20231011_150952'
px_ott = 0.00076 #mm per px
f0 = 0
f1 = 40
par_remapped = imgreg.load_registeredPar(tnpar)
par_filtered = th.comp_filtered_image(par_remapped,  d=px_ott, verbose=True, disp=False, freq2filter=(f0,f1))
par=par_filtered

par = np.ma.masked_array(par.data, mask)
dd = eye_ave.data -2*par.data


#### con maskera ad occhio di tigre !
auxmask = geo.draw_mask(np.ones(np.shape(eye_ave.data)),1030, 1000, 60,out=1)
mask2 =(-1* mask+1)*(-1*auxmask+1)
idaux = np.where(mask2==1)
mask[idaux]=1
dd1 = np.ma.masked_array(dd, mask)

figure(figsize=(10,5))
subplot(1,2,1); imshow(sum_mask*mask)
subplot(1,2,2); imshow(dd1); title('RMS = %.3e m' %(dd1.std()))

#### con maschera centrata nel max numero di img mediate con mascherina per togliere il centro
cc = [1250,1000,150]
nm = mask.copy()
nm = geo.draw_mask(np.ones(np.shape(eye_ave.data)),cc[0],cc[1],cc[2],out=1)
ndd = np.ma.masked_array(dd.data, nm)

figure(figsize=(10,5))
subplot(1,2,1); imshow(sum_mask*nm)
subplot(1,2,2); imshow(ndd); title('RMS = %.3e m' %(ndd.std()))


#### con maschera centrata in un punto dove ho la media fatta da un numero di immagini diverse (qui il raggio della maschera può essere a piacere)
cc = [1400, 1000, 150]
nm = mask.copy()
nm = geo.draw_mask(np.ones(np.shape(eye_ave.data)),cc[0],cc[1],cc[2],out=1)
ndd = np.ma.masked_array(dd.data, nm)

figure(figsize=(10,5))
subplot(1,2,1); imshow(sum_mask*nm)
subplot(1,2,2); imshow(ndd); title('RMS = %.3e m' %(ndd.std()))




show()

#analysis of the diff wrt PAR
img0 = ave.copy()
mmask = sum_mask.copy()
mmask = (mmask == np.max(mmask))
cir = geo.qpupil(mmask)
cir=[1400,1000]
rr = 300
mmask = (geo.draw_mask(np.zeros(shape(mmask)),cir[0],cir[1],rr,out=0))
imshow(mmask)
img = np.ma.masked_array(img0.data, -1*mmask+1)
img = th.removeZernike(img,[1,2,3])
imshow(img)

thepar = imgreg.image_remask(par_remapped, img)
thepar = th.removeZernike(thepar,[1,2,3])
imshow(thepar)
dd = img-2*thepar
dr= dd.std()
clf();imshow(dd, vmin=-vm,vmax=vm);title(dr)
rp = int(rr/2)
dd1 = dd[cir[0]-rp:cir[0]+rp,cir[1]-rp:cir[1]+rp]
img1 = img[cir[0]-rp:cir[0]+rp,cir[1]-rp:cir[1]+rp]
p1 = 2*thepar[cir[0]-rp:cir[0]+rp,cir[1]-rp:cir[1]+rp]

x, y = th.comp_psd(dd1, norm='ortho',verbose=True,d=px_ott)
plot(x, x*y,'x');xscale('log');yscale('log')
xi, yi = th.comp_psd(img1, norm='ortho',verbose=True,d=px_ott)
plot(xi, xi*yi,'o')
xp, yp = th.comp_psd(p1, norm='ortho',verbose=True,d=px_ott)
plot(xp, xp*yp)
legend(['Diff','PARinOTT','2Par'])


par_filtered = th.comp_filtered_image(par_remapped,  d=px_ott, verbose=True, disp=False, freq2filter=(f0,f1))
ww = img-2*par_filtered
ww=th.removeZernike(ww, [1,2,3])
imshow(ww)

