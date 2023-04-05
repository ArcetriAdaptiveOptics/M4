import sys
sys.path.append('/home/m4/towerbridge/lib/python3.8/site-packages/')
import plico_interferometer
import os
from astropy.io import fits as pyfits
from m4.ground.timestamp import Timestamp


#functions
def saveima(name, masked_image):
    pyfits.writeto(name, masked_image.data)
    pyfits.append(name, masked_image.mask.astype(int))

def acquire(N,bpath,x,y,ref=None):    
    fname = '{:04d}'.format(x)+'_{:04d}'.format(y)
    if ref is not None:
        fname=fname+ref
    base =os.path.join(bpath,fname)
    os.mkdir(base)
    for jj in range(N):
        print(jj)
        ima = interf.wavefront()
        saveima(os.path.join(base,'{:04d}'.format(jj))+'.fits',ima)


#code
interf = plico_interferometer.interferometer('192.168.22.79', 7300)


Bpath="/mnt/m4storage/Data/M4Data/OPTData/RefMirror/"
tn = Timestamp.now()

fold=os.path.join(Bpath,tn)
os.mkdir(fold)

x = '?'
y = '?'
ref = None
acquire(64,fold,x,y,ref)
