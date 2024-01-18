'''
Author: Luca Oggioni
The script imports and makes a comparison between the images of the RM reconstruct by AMOS and the one obtain on the OTT by us.
RM by AMOS -> stiching of 8 subapertures
RM by us -> average of many positions on the parabola  

'''

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
from m4.analyzers import requirement_analyzer as req_anal
import scipy
import samplerate

## import the data from fits file
filename="/home/m4/Desktop/immagini/RM_8pos.fits"

if thetype == 'fits':
    hduList = pyfits.open(filename)
    img = hduList[0].data
    
flatAMOS=np.ma.array(img)/2
flatAMOS.mask=np.isnan(img)
flatAMOS=th.removeZernike(flatAMOS,[1, 2, 3, 4, 5, 6])

filename="/home/m4/Desktop/immagini/FlatMirror.fits"
flatOTT=th.read_phasemap(filename)
flatOTT=np.rot90(flatOTT,2)
px_ott = 0.00076 #m/px
Dm=0.59; Rpx=int(Dm/px_ott/2)


flatOTT=flatOTT[440-390:440+390,440-390:440+390]

mm = geo.draw_mask(np.ones(np.shape(flatOTT.data)), np.shape(flatOTT)[0]/2, np.shape(flatOTT)[1]/2, Rpx, out=1)
flatOTT.mask[mm!=0]=1 

close("all")
opts = {'vmin': -5e-8, 'vmax': 5e-8, 'edgecolors': 'none'}
fig, (ax1, ax2) = subplots(1,2)
pc1=ax1.pcolormesh(flatAMOS, **opts); fig.colorbar(pc1); ax1.set_aspect('equal','box')
pc2=ax2.pcolormesh(flatOTT, **opts); fig.colorbar(pc2); ax2.set_aspect('equal','box')
ax1.set_title('AMOS, RMS: {:.3}'.format(flatAMOS.std()))
ax2.set_title('OTT, RMS: {:.3}'.format(flatOTT.std()))
fig.tight_layout()
show()

###### remapping and subtraction
# find the best angle of rotation
x,y=np.shape(flatOTT)
rr=[]
ang=np.arange(10,20,0.1)
for t in ang:
    tmp=scipy.ndimage.rotate(flatOTT, angle=t)
    x1,y1=np.shape(tmp)
    tmp=tmp[int((y1-y)/2):int((y1+y)/2),int((y1-y)/2):int((y1+y)/2)]
    
    #figure();imshow(tmp); show()
    
    flatOTT2 = scipy.ndimage.zoom(tmp.data, np.shape(flatAMOS)[0]/np.shape(flatOTT)[0])  
    mm = scipy.ndimage.zoom(flatOTT.mask, np.shape(flatAMOS)[0]/np.shape(flatOTT)[0])  
    flatOTT2=np.ma.array(flatOTT2); flatOTT2.mask=flatAMOS.mask
    
    res=flatOTT2-flatAMOS
    
    rr.append(np.std(res))
    
figure(); plot(ang,rr); show()

##rotate and subtract 
t=13.4
tmp=scipy.ndimage.rotate(flatOTT, angle=t)
x1,y1=np.shape(tmp)
tmp=tmp[int((y1-y)/2):int((y1+y)/2),int((y1-y)/2):int((y1+y)/2)]

#figure();imshow(tmp); show()

flatOTT2 = scipy.ndimage.zoom(tmp.data, np.shape(flatAMOS)[0]/np.shape(flatOTT)[0])  
mm = scipy.ndimage.zoom(flatOTT.mask, np.shape(flatAMOS)[0]/np.shape(flatOTT)[0])  
flatOTT2=np.ma.array(flatOTT2); flatOTT2.mask=flatAMOS.mask

res=flatOTT2-flatAMOS

close("all")
opts = {'vmin': -5e-8, 'vmax': 5e-8, 'edgecolors': 'none'}
fig, (ax1, ax2) = subplots(1,2)
pc1=ax1.pcolormesh(flatAMOS, **opts); fig.colorbar(pc1); ax1.set_aspect('equal','box')
pc2=ax2.pcolormesh(flatOTT2, **opts); fig.colorbar(pc2); ax2.set_aspect('equal','box')
ax1.set_title('AMOS, RMS: {:.3}'.format(flatAMOS.std()))
ax2.set_title('OTT, RMS: {:.3}'.format(flatOTT2.std()))
fig.tight_layout()
show()

figure(); imshow(res); colorbar(); clim(-5e-8,5e-8);title('res, RMS: {:.3}'.format(res.std()))
show()

## check 

im_flat=flatAMOS

px_ott = 0.001*532/600 #mm per px

l=int(0.05/px_ott) #semilato campionamento 
step=int(0.01/px_ott) #step di campionamento

start_px=100
stop_px=400


cutx=[]
cuty=[]
cutxy=[]
for i in np.arange(start_px,stop_px,step):
    plot(266,i,'x')
    cutx.append(im_flat[266-l:266+l,i:i+l*2])
    cuty.append(im_flat[i:i+l*2,266-l:266+l])
    cutxy.append(im_flat[266-l:266+l,i:i+l*2])
    cutxy.append(im_flat[i:i+l*2,266-l:266+l])
    


cutx=np.ma.masked_array(cutx)
patchx=np.ma.mean(cutx,axis=0)
cuty=np.ma.masked_array(cuty)
patchy=np.ma.mean(cuty,axis=0)
cutxy=np.ma.masked_array(cutxy)
patchxy=np.ma.mean(cutxy,axis=0)

opts = {'vmin': -0.8e-8, 'vmax': 0.8e-8, 'edgecolors': 'none'}
fig, (ax1, ax2, ax3) = subplots(1,3)
pc1=ax1.pcolormesh(patchx, **opts); fig.colorbar(pc1); ax1.set_aspect('equal','box')
pc2=ax2.pcolormesh(patchy, **opts); fig.colorbar(pc2); ax2.set_aspect('equal','box')
pc3=ax3.pcolormesh(patchxy, **opts); fig.colorbar(pc3); ax3.set_aspect('equal','box')
ax1.set_title('shift along x, RMS: {:.3}'.format(patchx.std()))
ax2.set_title('shift along y, RMS: {:.3}'.format(patchy.std()))
ax3.set_title('shift along x and y, RMS: {:.3}'.format(patchxy.std()))
fig.tight_layout()
show()
# figure(); imshow(patchx); show()