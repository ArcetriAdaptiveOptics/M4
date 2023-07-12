import sys
sys.path.append('/home/m4/towerbridge/lib/python3.8/site-packages/')
import plico_interferometer
import os
import time
from astropy.io import fits as pyfits
from m4.ground.timestamp import Timestamp


#functions
def saveima(name, masked_image):
    pyfits.writeto(name, masked_image.data)
    pyfits.append(name, masked_image.mask.astype(int))

def acquire(N,bpath,x,y,ref=None,delay=None):    
    fname = '{:04d}'.format(x)+'_{:04d}'.format(y)
    if ref is not None:
        fname=fname+ref
    base =os.path.join(bpath,fname)
    os.mkdir(base)
    for jj in range(N):
        print(jj)
        if delay!= None:
            print("waiting "+str(delay)+"s ...")
            time.pause(delay)
            
        ima = interf.wavefront()
        saveima(os.path.join(base,'{:04d}'.format(jj))+'.fits',ima)
        

# ['',
#  '/home/m4/anaconda3/lib/python38.zip',
#  '/home/m4/anaconda3/lib/python3.8',
#  '/home/m4/anaconda3/lib/python3.8/lib-dynload',
#  '/home/m4/anaconda3/lib/python3.8/site-packages',
#  '/home/m4/anaconda3/lib/python3.8/site-packages/IPython/extensions',
#  '/home/m4/towerbridge/lib/python3.8/site-packages/',
#  '/home/m4/git/M4']

#code
interf = plico_interferometer.interferometer('192.168.22.79', 7300)


Bpath="/mnt/m4storage/Data/M4Data/OPTData/RefMirror/"
tn = Timestamp.now()

fold=os.path.join(Bpath,tn)
os.mkdir(fold)

# x = '?'
# y = '?'
# ref = None
# acquire(10,fold,x,y,ref)
