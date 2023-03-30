##
# script per analizzare le immagini full frame e calcolare il gradiente x,y da confrontare con le misure di AMOS

from astropy.io import fits as pyfits
from m4.mini_OTT import timehistory as th
from matplotlib import pyplot as plt
from m4.ground import zernike
import numpy as np
from m4 import noise
from m4.configuration import config_folder_names as fold_name
from m4.ground import rebinner
import os

# definisco la maschera ma
nx=425*4
ny=425*4
x = np.linspace(1, nx, nx)
y = np.linspace(1, ny, ny)
xv, yv = np.meshgrid(x, y)
x0=215.5*4;y0=215.5*4;R=(26.5+8)*4;
mask1=(xv-x0)**2+(yv-y0)**2<R**2
x0=215.5*4;y0=214.5*4;R=(196.5-10)*4;
mask2=(xv-x0)**2+(yv-y0)**2>R**2
ma=mask1 | mask2

##
# calcolo l'average, tolgo zernike, 
# plotto rms della differenza tra i singoli frame e la media
# plotto l'immagine dell'average 

tn=[]
#tn.append('20230303_095811')
tn.append('20230303_102203')

idx=[0, 99]

fl = th.fileList(tn[0])
q0=th.averageFrames(idx[0],idx[1],fl)
#q0.mask=ma  
coeff, mat = zernike.zernikeFit(q0, [5,6])
zernCoef=coeff
q1=th.removeZernike(q0,[1,2,3,4,7,8])

rr = []
for i in range(idx[0],idx[1]):
    qi = th.removeZernike(th.frame(i,fl)-q0,[1,2,3,4,7,8])
    rr.append(qi.std())



plt.figure(1); plt.clf()
plt.plot(rr, '-o'); plt.title('RMS-ref'); plt.legend(tn)

plt.figure(2); plt.clf()
plt.imshow(q1)
plt.title(tn[0]+' RMS={:.2e}'.format(np.std(q1))+'\n z5={:.2e}'.format(zernCoef[0])+', z6={:.2e}'.format(zernCoef[1]));  
plt.colorbar()    

plt.show()

#
# calcolo la differenza shiftando l'immagine di N pixel
# rebinna la differenza

#plt.close('all')

dpix = 1
dd = q1-np.roll(q1,(dpix,dpix),axis=(0,1))

plt.figure(3);plt.clf();plt.imshow(dd);plt.colorbar();
plt.title('PAR gradient \n'+tn[0]+' RMS={:.2e}'.format(np.std(dd))+' shift {:d} pixel'.format(dpix)); 
plt.clim(-8e-9,8e-9)
plt.show();

K=2
new_shape=(dd.shape[0]/K, dd.shape[1]/K )
dd_b=rebinner.rebin2DArray(dd, new_shape, sample=False)

plt.figure(4);plt.clf();plt.imshow(dd_b);plt.colorbar();
plt.title('PAR gradient_rebinned, 2pix \n'+tn[0]+' RMS={:.2e}'.format(np.std(dd_b))); 
plt.clim(-8e-9,8e-9)
plt.show();

fits_file_name = os.path.join(dove, name)
            pyfits.writeto(fits_file_name, masked_ima.data)
            pyfits.append(fits_file_name, masked_ima.mask.astype(int))
