#script for the comparison of the PAR shape before-after coating
from m4.ground import zernike as zern
from m4.ground import geo as geo
from matplotlib.pyplot import *
from astropy.io import fits
from m4.ground import timestamp
from m4.ground import read_data
from m4.mini_OTT import timehistory as th
import numpy as np
import os



#finding the dataset...
img2flat5 = ['20221216_104500-16ave.4D','20221216_104000-16ave.4D' ]  # corresponda to p00bis
img2flat4bis = ['20220804_155900.4D','20220804_160000.4D','20220804_160100.4D','20220804_160200.4D','20220804_160300.4D']  # no cmd applied, for repatab./stability vs imgflat4

def read_frame4d(tn, name):
    base = '/mnt/m4storage/Data/M4Data/OPTData/PARTest/'
    name = base+tn+'/'+name
    img=read_data.read_phasemap(name)
    return img



tn0 = '20220802'
name0 = ['20220802_123600.4D','20220802_123700.4D','20220802_123800.4D','20220802_123900.4D']
name0 = ['20220804_155900.4D','20220804_160000.4D','20220804_160100.4D','20220804_160200.4D','20220804_160300.4D']
tn1 = '20221219'
name1 = ['20221219_101200.4D','20221219_101300.4D']

img0 = read_frame4d(tn0,name0[0])
img0 = removeZernike(img0,[1,2,3,4,5,6,7,8])
clf();imshow(img0);colorbar();title('PAR before coating')
img0.std()

m0 = img0.mask

img1 = read_frame4d(tn1,name1[0])
m1 = (np.invert(img1.mask))
cir = geo.qpupil(m1)
m1 = geo.draw_mask(m1, cir[0],cir[1],cir[2])
inrad = 500
m2 = geo.draw_mask(m1, cir[0],cir[1],inrad,out=1)

img1 = np.ma.masked_array(img1, np.invert(m2))
img1 = removeZernike(img1,[1,2,3,4,5,6,7,8])
clf()
imshow(img1)
img1.std()

tns = '20230127_163227'
tns = '20230223_224400'
fl = th.fileList(tns)
img2 = th.averageFrames(0,100,fl)
cc,cmat = th.zernike.zernikeFit(img2,[1,2,3,4,5,6,7,8,9,10])
img2 = th.removeZernike(img2, [1,2,3,4,7,8])
img2.std()
img2 = th.removeZernike(img2, [1,2,3,4,5,6,7,8])
img2.std()*2
clf()
imshow(img2)
clf();imshow(img2);colorbar();title('PAR after coating')

dpix = 1
dd = img2-np.roll(img2,(dpix,dpix),axis=(0,1))
clf();imshow(dd);colorbar();title('PAR gradient, 2pix');
dd.std()*2

#ritagliare sulla CA


