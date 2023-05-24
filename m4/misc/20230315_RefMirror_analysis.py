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
from m4 import noise

base = '/mnt/m4storage/Data/M4Data/OPTData/RefMirror/'
#tn   = '20230331_124514'
#tn   = '20230418_100251'
tn = '20230418_113032'
avefold=os.path.join(base,tn+'_Averages')
if not os.path.exists(avefold):
    os.mkdir(avefold)

### stich analysis: produce averages

listfolder=os.listdir(base+tn)

qq=1
for fold in listfolder:
    print(qq); qq=qq+1
    plt.close('all')  
    imglist = os.listdir(base+tn+'/'+fold)
    
    st = []
    somma = []
    for i in imglist:
        img=th.read_phasemap(base+tn+'/'+fold+'/'+i)
        maschera=(np.invert(img.mask));
        somma.append(np.sum(maschera));
    #    plt.clf(); plt.imshow(img); plt.pause(0.5); plt.show()
    
        
    #plt.figure();plt.plot(somma); plt.show()
    
    c= np.median(somma)
    idgood = somma > c*0.9
    imglist=list(compress(imglist, idgood))
    for i in imglist:
        img=th.read_phasemap(base+tn+'/'+fold+'/'+i)
        img = th.removeZernike(img,[1,2,3])
        st.append(img.std())
    #plt.clf(); plt.imshow(img); plt.pause(0.5); plt.show()
    
    #plt.figure();plt.plot(st); plt.show()   
    
    c= min(st)
    idgood = np.where(st < c*3)
    idgood = st < c*3
    imglist=list(compress(imglist, idgood))       
    
    av=aveframe3(os.path.join(base,tn,fold),imglist)
    plt.figure(); plt.imshow(av);  colorbar(); plt.show();    
        
    save_phasemap(avefold, fold, av)



# ## guarda un'immagine alla volta
# i=0
# ##
# plt.close('all')
# plt.figure()
# img = th.read_phasemap(base+tn+'/'+fold+'/'+list[i])
# img= th.removeZernike(img,[1,2,3])
# print(img.std())
# plt.imshow(img)
# plt.show(); 
# i=i+1
# #
# # movie(tn,fold)


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

def aveframe3(fold,imglist):
    
    cou = 0
    imgcube=[]
    for ii in imglist:
        fname = os.path.join(fold,ii)
        img = th.read_phasemap(fname)
        imgcube.append(img)
    imgAve=np.ma.mean(imgcube, axis=0)
    return imgAve

def fit_circle_2d(x, y, w=[]):
    
    A = np.array([x, y, np.ones(len(x))]).T
    b = x**2 + y**2
    
    # Modify A,b for weighted least squares
    if len(w) == len(x):
        W = diag(w)
        A = dot(W,A)
        b = dot(W,b)
    
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

def save_phasemap(location, file_name, masked_image):
    """
    Parameters
    ----------
    location: string
        measurement file path
    file_name: string
        measuremnet fits file name
    masked_image: numpy masked array
        data to save
    """
    fits_file_name = os.path.join(location, file_name)
    pyfits.writeto(fits_file_name, masked_image.data)
    pyfits.append(fits_file_name, masked_image.mask.astype(int))



##
listfolder=os.listdir(base+tn)
listfolder.sort()

tau_vector = np.arange(1,5,1)
path_series = '/mnt/m4storage/Data/M4Data/OPTData/RefMirror/'

dfpath=os.path.join(path_series,tn,listfolder[1])

noise.convection_noise(dfpath, tau_vector)