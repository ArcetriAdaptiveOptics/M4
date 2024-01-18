import numpy as np
from m4.misc import par_registration_lib as parreg
from m4.ground import geo
from m4.mini_OTT import timehistory as th
from importlib import reload

#init
#PAR meas with CGH
cgh_tn_marker = '20230428_093328'; off_cgh_marker = [180,180]
cgh_tn_img = '20230428_164257'; off_cgh_img = [180,180]
#OTT meas with RM
ott_tn_marker = '20230707_142156'; off_ott_marker = [24,0]#was 0,24
ott_tn_img = '20230715_113735'; off_ott_img = [600,600]

#CGH vs OTT markers match, at Center_view
mark_cgh_list = np.array([8,12,11,6,13,18,17,14,7])
mark_ott_list= np.array([2,4,3,1,6,8,7,5,0])

#definitio of image dimensions
ott_img_diam = 0.6/2
cgh_img_diam = 1.44/2
ottdpix = 800
zlist = np.arange(10)+1
fl = th.fileList(tn)
n=400
img0=th.averageFrames(0,n,fl)
img1=th.averageFrames(n+1,n+n+1,fl)
dd = img1-img0
dd = th.removeZernike(dd,zlist)
imshow(dd, vmin=-3e-9,vmax=3e-9)

w = quick243(dd,pixs)
print(mean(w))

cir=geo.qpupil(-1*img.mask+1)
pixs = cir[2]/ott_img_diam

cc = th.cubeFromList(fl[0:200:2])
cc1 = []
for i in cc:
    g=th.removeZernike(i,zlist)
    cc1.append(g)
cc1 = np.ma.masked_array(cc1)

def quick243(img, pixs):  #pixs = [pix/m]
    img1 = crop_frame(img)
    dd = 0.03
    pp =int((pixs * dd))#.astype(int)
    ss = np.array(img1.shape)
    ss = (ss/pp).astype(int)
    ww = np.zeros(ss)
    for ii in np.arange (0,ss[0]):
        for jj in np.arange(0,ss[1]):
            kk = img1[ii*pp:ii*pp+pp,jj*pp:jj*pp+pp]
            st = kk.std()
            if st != np.nan:
                ww[ii,jj] = st
    mask = (ww == 0)
    ww = np.ma.masked_array(ww.data,mask)
    return ww


def crop_frame(imgin):
    cir = geo.qpupil(-1*imgin.mask+1)
    cir = np.array(cir[0:3]).astype(int)
    img = imgin.data[cir[0]-cir[2]:cir[0]+cir[2]+1,cir[1]-cir[2]:cir[1]+cir[2]+1]
    m = imgin.mask[cir[0]-cir[2]:cir[0]+cir[2]+1,cir[1]-cir[2]:cir[1]+cir[2]+1]
    img = np.ma.masked_array(img, m)
    return img
