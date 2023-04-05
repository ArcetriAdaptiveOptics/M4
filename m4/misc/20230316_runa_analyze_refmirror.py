from m4.ground import zernike as zern
from m4.ground import geo as geo
from matplotlib.pyplot import *
from astropy.io import fits as pyfits
from m4.ground import timestamp
from m4.ground import read_data
from m4.mini_OTT import timehistory as th
from itertools import compress
import numpy as np
from m4.ground import zernike
from matplotlib import pyplot as plt
import os


### functions
def gimmethelist(tn, fold):
    flist = os.listdir(base+tn+'/'+fold+'/')
    nf = len(flist)
    f2 = []
    for i in flist:
        f2.append(base+tn+'/'+fold+'/'+i)
    return f2

def readframe(name):
    img= th.read_phasemap(name)
    return img

def movie(tn,fold,time=0.5):
    imglist = os.listdir(base+tn+'/'+fold)
    nf = len(imglist)
    plt.figure()
    for i in range(nf):
        fname = base+tn+'/'+fold+'/'+imglist[i]
        img = th.read_phasemap(fname)  
        plt.clf()
        plt.imshow(img)
        plt.axis('off')
        #plt.colorbar()
        plt.pause(time)
        plt.show()

def aveframe(tn,fold, thr = 0.9):
    imglist = os.listdir(base+tn+'/'+fold)

    nf = len(imglist)
    cou = 0
    for i in range(nf):
        fname = base+tn+'/'+fold+'/'+imglist[i]
        img = th.read_phasemap(fname)
        ss = np.shape(img) 
        if i == 0:
            imgmaster = np.ma.masked_array(np.zeros(ss),np.ma.make_mask_none(ss)) 
        cir = geo.qpupil(img.mask)
        npp = np.pi*cir[2]**2
        npv = np.size(img.mask)-np.sum(img.mask.astype(int))
        isgood = npv > npp*thr
        if isgood == True:
            print(i)
            cou = cou+1
            imgmaster = np.ma.masked_array(img.data+imgmaster.data, np.ma.mask_or(imgmaster.mask, img.mask))
    imgmaster = np.ma.masked_array(imgmaster.data/cou,imgmaster.mask)
    return imgmaster
    

def aveframe2(tn,fold, idgood=np.array([])):
    
    imglist = os.listdir(base+tn+'/'+fold)
    nf = len(imglist)
    if idgood.any()==False:
        idgood=np.full((nf, 1), True)

    
    cou = 0
    imgcube=[]
    for i in np.arange(nf):
        if idgood[i]==True:
            fname = base+tn+'/'+fold+'/'+imglist[i]
            img = th.read_phasemap(fname)
            imgcube.append(img)
    imgAve=np.ma.mean(imgcube, axis=0)
    return imgAve


def fit_circle_2d(x, y, w=[]):
    
    A = np.array([x, y, np.ones(len(x))]).T
    b = x**2 + y**2
    
    # Modify A,b for weighted least squares
    if len(w) == len(x):
        W = np.diag(w)
        A = np.dot(W,A)
        b = np.dot(W,b)
    
    # Solve by method of least squares
    c = np.linalg.lstsq(A,b,rcond=None)[0]
    
    # Get circle parameters from solution c
    xc = c[0]/2
    yc = c[1]/2
    r = np.sqrt(c[2] + xc**2 + yc**2)
    return xc, yc, r

def std_check(flist, thr):
    nf = len(flist)
    st = []
    for i in flist:
        img=th.read_phasemap(i)
        img = th.removeZernike(img,[1,2,3])
        st.append(img.std())
    c= min(st)
    idgood = np.where(st < c*thr)
    return idgood, st


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

#code
base = '/mnt/m4storage/Data/M4Data/OPTData/RefMirror/'
tn   = '20230329_091725'
fdlist = os.listdir(base+tn)
reffold = '0000_0000_ref0'
reflist = os.listdir(base+tn+'/'+reffold)
f=0
imglist = os.listdir(base+tn+'/'+fdlist[f])
i=0
fname = base+tn+'/'+fdlist[0]+'/'+imglist[i]

refname = base+tn+'/'+fdlist[0]+'/'+imglist[1]
ref = th.read_phasemap(refname)

## reference analysis
plt.close('all')

reffold = '0000_0000_ref0'
reflist = os.listdir(base+tn+'/'+reffold)

nf = len(reflist)
st = []
somma = []
for i in reflist:
    img=th.read_phasemap(base+tn+'/'+reffold+'/'+i)
    maschera=(np.invert(img.mask));
    somma.append(np.sum(maschera));
#    plt.clf(); plt.imshow(img); plt.pause(0.5); plt.show()

    
plt.figure();plt.plot(somma); plt.show()

c= np.median(somma)
idgood = somma > c*0.9
reflist=list(compress(reflist, idgood))
plt.figure()
for i in reflist:
    img=th.read_phasemap(base+tn+'/'+reffold+'/'+i)
    img = th.removeZernike(img,[1,2,3])
    st.append(img.std())
#    plt.clf(); plt.imshow(img); plt.pause(0.5); plt.show()

plt.figure();plt.plot(st); plt.show()   

c= min(st)
idgood = np.where(st < c*3)
idgood = st < c*3
reflist=list(compress(reflist, idgood))

av=aveframe2(tn, reffold)
plt.figure(); plt.imshow(av);  colorbar(); plt.show();    
        


    
##
movie(tn,reffold)
#av=aveframe(tn, reffold, 0.5)
av=aveframe2(tn, reffold)
plt.imshow(av); plt.show()
##
# definisco la maschera ma

nx=990
ny=998
# maskk=np.invert(av.mask).astype(int)
# cir = geo.qpupil(np.invert(av.mask).astype(int))
# x0=cir[0]
# y0=cir[1]
# r=cir[2]
# xx=cir[3]
# yy=cir[4]
x = np.linspace(1, nx, nx)
y = np.linspace(1, ny, ny)
xx, yy = np.meshgrid(x, y)
x0=520;y0=520;R=173;
mask1=(xx-x0)**2+(yy-y0)**2>R**2
ma=mask1 
#
av.mask=ma
coeff, mat = zernike.zernikeFit(av, [1,2,3,4,5,6,7,8])
print(coeff)
av_zern=th.removeZernike(av,[1,2,3,4])
plt.figure(); plt.imshow(av_zern); plt.colorbar(); plt.show()

### stich analysis
plt.close('all')

fold = '0002_0200'
imglist = os.listdir(base+tn+'/'+fold)

nf = len(imglist)
st = []
for i in imglist:
    img=th.read_phasemap(base+tn+'/'+fold+'/'+i)
    img = th.removeZernike(img,[1,2,3])
    st.append(img.std())

plt.figure(); plt.plot(st); plt.show()

c= min(st)
idgood = np.where(st < c*3)
idgood = st < c*3

av=aveframe2(tn, fold, idgood)
plt.figure(); plt.imshow(av);  colorbar(); plt.show();


i=0
##
plt.close('all')
plt.figure()
img = th.read_phasemap(base+tn+'/'+fold+'/'+list[i])
img= th.removeZernike(img,[1,2,3])
print(img.std())
plt.imshow(img)
plt.show(); 
i=i+1
#
# movie(tn,fold)

## definisco la maschera ma
# p1=plt.ginput(3)
# xp=np.array([p1[0][0], p1[1][0], p1[2][0]])
# yp=np.array([p1[0][1], p1[1][1], p1[2][1]])
# [xc1, yc1, r1]=fit_circle_2d(xp,yp)
# 
# p1=plt.ginput(3)
# xp=np.array([p1[0][0], p1[1][0], p1[2][0]])
# yp=np.array([p1[0][1], p1[1][1], p1[2][1]])
# [xc2, yc2, r2]=fit_circle_2d(xp,yp)
# 
# nx=990
# ny=998
# x = np.linspace(1, nx, nx)
# y = np.linspace(1, ny, ny)
# xv, yv = np.meshgrid(x, y)
# mask1=(xv-xc1)**2+(yv-yc1)**2>r1**2
# mask2=(xv-xc2)**2+(yv-yc2)**2>r2**2
# 
# ma=mask1 | mask2
# av.mask=ma | av.mask
# ##

coeff, mat = zernike.zernikeFit(av, [1,2,3,4,5,6,7,8])
print(coeff)
av_zern=th.removeZernike(av,[1,2,3,4])
plt.figure(); plt.imshow(av_zern); plt.colorbar(); plt.show()
