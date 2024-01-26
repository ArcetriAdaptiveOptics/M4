import numpy as np
from importlib import reload
import os
from astropy.io import fits as pyfits
from m4.ground import zernike as zern
from m4.mini_OTT import timehistory as th
from m4.ground import geo, rebinner



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

base = '/mnt/m4storage/Data/M4Data/OPTData/RefMirror/'
tn = '20230315_155719'

stthresh = 2
maskthresh=0.9


def rebin(img0,factor=8):
    img = img0.data
    mask = img0.mask
    ss = np.shape(img)
    aa = (int(ss[0]/factor)*factor, int(ss[1]/factor)*factor)
    imgcrop = img[0:aa[0],0:aa[1]]
    maskcrop = mask[0:aa[0],0:aa[1]]
    ss1 = (int(aa[0]/factor), int(aa[1]/factor))
    img1 =  rebinner.rebin2DArray(imgcrop,ss1,sample=False)
    mask1 = rebinner.rebin2DArray(maskcrop,ss1,sample=False) 
    img = np.ma.masked_array(img1,mask1)
    return img

    

def fold2pos(tn):
    foldlist = sorted(os.listdir(base+tn))
    p0 = []
    p1 = []
    fl = []
    for i in foldlist:
        ch = i.split('_')
        if len(ch) == 2:
            p0.append(int(ch[0])) 
            p1.append(int(ch[1]))       
            fl.append(i)
    npos = len(p0)
    ppos = np.zeros([npos,2])
    ppos[:,0] = p0
    ppos[:,1] = p1
    return ppos, fl

def reffold(tn):
    foldlist = os.listdir(base+tn)
    reflist = []
    for i in foldlist:
        ch = i.split('_')
        if len(ch) == 3:
            reflist.append(i)
    return reflist

def getrefdata(tn, refdiam=.1):
    flref = reffold(tn)
    img = patchdata(tn, flref[0],thrsurf=2, thrpix=0.8)
    cavity = th.removeZernike(img,[1,2,3])
    cir=geo.qpupil(np.invert(img.mask))
    pixscale = cir[2]*2/refdiam
    return pixscale





def gimmethelist(tn, fold):
    flist = os.listdir(base+tn+'/'+fold+'/')
    nf = len(flist)
    f2 = []
    for i in flist:
        f2.append(base+tn+'/'+fold+'/'+i)
    return f2


def getcube(flist, rebfactor):
    nf = len(flist)
    imgcube = []
    for i in flist:
        img = th.read_phasemap(i)
        img = rebin(img,rebfactor)
        imgcube.append(img)
    return imgcube
    
def std_check(thecube, thr, removeZern=0):
    nc = len(thecube)
    st = []
    for i in thecube:
        if removeZern == 1:    
            img = th.removeZernike(i,[1,2,3])
        else:
            img = i
        st.append(np.ma.std(img))
    vv = np.min(st)
    idgood = []
    for i in range(nc):
        if st[i] < vv*thr:
            idgood.append(thecube[i])
    #idgood = [x for x in range(nc) if x < vv*thr]
    #idgood = np.where(st > vv*thr)
    return idgood

def mask_check(thecube, thr):
    nc = len(thecube)
    st = []
    for i in thecube:
        mi = np.invert(i.mask).astype(int)
        st.append(np.sum(mi))
    vv = np.max(st)
    idgood = []
    for i in range(nc):
        if st[i] > vv*thr:
            idgood.append(thecube[i])
    #idgood = [x for x in range(nc) if x > vv*thr]
    #idgood = np.where(st < vv*thr)
    return idgood

def patchdata(tn,fold, thrsurf=2, thrpix=0.8, rebinfactor=8):   
    fl = gimmethelist(tn, fold)
    cc = getcube(fl,rebinfactor)
    cc = mask_check(cc, thrpix)
    print(len(cc))
    cc = std_check(cc, thrsurf)
    print(len(cc))
    img=np.ma.mean(cc, axis=0)
    return img



def tndata(tn, thrsurf=2, thrpix=0.8, rebinfactor=8):
    ppos,plist = fold2pos(tn)
    nf = len(plist)
    patch = []
    maskvec = []
    for i in plist:
        print(i)
        pp = patchdata(tn, i, thrsurf, thrpix, rebinfactor=8)
        mp = np.invert(pp.mask).astype(int)
        patch.append(pp)
        maskvec.append(mp)
    return patch, maskvec, ppos


def tndataref(tn, thrsurf=2, thrpix=0.8, rebinfactor=8):
    plist = reffold(tn)
    nf = len(plist)
    patch = []
    maskvec = []
    for i in plist:
        print(i)
        pp = patchdata(tn, i, thrsurf, thrpix, rebinfactor=8)
        mp = np.invert(pp.mask).astype(int)
        patch.append(pp)
        maskvec.append(mp)
    return patch, maskvec

def tndataall(tn, thrsurf=2, thrpix=0.8, rebinfactor=8):  
    foldlist = os.listdir(base+tn)
    foldlist = np.sort(foldlist)
    reflist = []
    imglist = []
    for i in foldlist:
        ch = i.split('_')
        ch = ch[-1]    
        if ch[0:2] == 're':
            reflist.append(i)
        if ch[0:2] == 'im':
            imglist.append(i)
    refpatch = []
    refmaskvec = []
    for i in reflist:
        print(i)
        pp = patchdata(tn, i, thrsurf, thrpix, rebinfactor=8)
        mp = np.invert(pp.mask).astype(int)
        refpatch.append(pp)
        refmaskvec.append(mp)        
    imgpatch = []
    imgmaskvec = []
    for i in imglist:
        print(i)
        pp = patchdata(tn, i, thrsurf, thrpix, rebinfactor=8)
        mp = np.invert(pp.mask).astype(int)
        imgpatch.append(pp)
        imgmaskvec.append(mp)
    #imgpatch = np.ma.masked_array(imgpatch)
    #refpatch = np.ma.masked_array(refpatch)

    return imgpatch, imgmaskvec, refpatch, refmaskvec



def erodemask(maskvec,imgvec,dpix=1):
    for ii in range(len(imgvec)):
        maskvec[ii] = maskvec[ii] * np.roll(maskvec[ii],(dpix,dpix),axis=(0,1))
        maskvec[ii] = maskvec[ii] * np.roll(maskvec[ii],(-1*dpix,-1*dpix),axis=(0,1))
        imgvec[ii].mask = -1*maskvec[ii]+1
    return imgvec, maskvec



def refmirr_fullmask(maskvec, pos, pixscale):
    #pos shall be nposx2 in [m]
    #pixscale shall be [pix/m]
    #shall be run on masks only
    npos = len(pos)
    pixpos = (pos*pixscale).astype(int)
    imgsize = np.array(np.shape(maskvec[0])).astype(int)
    imgfsize=imgsize+pixpos.max(0)
    imgf=np.ones([int(imgfsize[0]),int(imgfsize[1])])
    for i in range(len(pixpos)):
        #imgf[pixpos[i,0]:pixpos[i,0]+imgsize[0],pixpos[i,1]:pixpos[i,1]+imgsize[1]] = imgf[pixpos[i,0]:pixpos[i,0]+imgsize[0],pixpos[i,1]:pixpos[i,1]+imgsize[1]] +maskvec[:,:,i]
        imgf[pixpos[i,0]:pixpos[i,0]+imgsize[0],pixpos[i,1]:pixpos[i,1]+imgsize[1]] = imgf[pixpos[i,0]:pixpos[i,0]+imgsize[0],pixpos[i,1]:pixpos[i,1]+imgsize[1]] *(-1*maskvec[i]+1)
    
    imgf = -1*imgf+1
    return imgf

def refmirr_fullframes(maskvec,imgvec, pos, pixscale):
    #pos shall be nposx2 in [m]
    #pixscale shall be [pix/m]
    #shall be run on masks only
    npos = len(pos)
    pixpos = (pos*pixscale).astype(int)
    imgsize = np.array(np.shape(maskvec[0])).astype(int)
    imgfsize=imgsize+pixpos.max(0)
    maskf = np.zeros([int(imgfsize[0]),int(imgfsize[1])])
    maskf = np.ones([int(imgfsize[0]),int(imgfsize[1])])
    imgf=np.zeros([int(imgfsize[0]),int(imgfsize[1])])
    for i in range(len(pixpos)):
        #imgf[pixpos[i,0]:pixpos[i,0]+imgsize[0],pixpos[i,1]:pixpos[i,1]+imgsize[1]] = imgf[pixpos[i,0]:pixpos[i,0]+imgsize[0],pixpos[i,1]:pixpos[i,1]+imgsize[1]] +maskvec[:,:,i]
        maskf[pixpos[i,0]:pixpos[i,0]+imgsize[0],pixpos[i,1]:pixpos[i,1]+imgsize[1]] = maskf[pixpos[i,0]:pixpos[i,0]+imgsize[0],pixpos[i,1]:pixpos[i,1]+imgsize[1]] *(-1*maskvec[i]+1)# +maskvec[i]
    maskf = -1*maskf +1
    imgvecout = []
    for i in range(len(pixpos)):
        imgf=np.ma.zeros([int(imgfsize[0]),int(imgfsize[1])])
        imgf[pixpos[i,0]:pixpos[i,0]+imgsize[0],pixpos[i,1]:pixpos[i,1]+imgsize[1]] = imgvec[i]
        imgf = imgf.data
        maskf = np.ones([int(imgfsize[0]),int(imgfsize[1])])
        maskf[pixpos[i,0]:pixpos[i,0]+imgsize[0],pixpos[i,1]:pixpos[i,1]+imgsize[1]] = maskf[pixpos[i,0]:pixpos[i,0]+imgsize[0],pixpos[i,1]:pixpos[i,1]+imgsize[1]] *(-1*maskvec[i]+1)
        imgf = np.ma.masked_array(imgf,maskf)
        imgvecout.append(imgf)
    return imgvecout

def merge(img1,img2):
    mask1 = -np.fix(img1.mask)+1
    mask2 = -np.fix(img2.mask)+1
    img1 = img1.data*mask1
    img2 = img2.data*mask2
    weight_mask = mask1 + mask2
    #idx = np.where(weight_mask != 0)
    master_mask = weight_mask*0.
    #master_mask[idx] = 1
    master_mask[weight_mask != 0] = 1

    master_img = np.ma.masked_array((img1 + img2),mask=master_mask==0)
    weight_img = (master_img/weight_mask)

    return weight_img

def coalign(img1, img2, dmmask, zlist2fit=[1,2,3,4], zlist2adjust=[1,2,3,4]):
    # img1 e img2 sono immagini grandi (img2dm)
    zlist2fit = np.array(zlist2fit)
    zlist2adjust = np.array(zlist2adjust)
    mask1 = -1*img1.mask +1
    mask2 = -1*img2.mask +1
    mask_inter = mask1*mask2
    img11 = np.ma.masked_array(img1*mask_inter, -1* mask_inter+1)
    img22 = np.ma.masked_array(img2*mask_inter, -1*mask_inter+1)

    cc1, zmat1 = zern.zernikeFitAuxmask(img11, dmmask, zlist2fit)
    cc2, zmat2 = zern.zernikeFitAuxmask(img22, dmmask, zlist2fit)
#    print('Coeff 0'); print(cc1)
#    print('Coeff 1'); print(cc1)

    dcc = cc2 - cc1
    dcc2adjust = dcc[zlist2adjust - 1]
#    dcc2adjust[0]=0
    tip2add = zernike_surf(zlist2adjust, -1*dmmask+1,dcc2adjust)*mask2
    img222 = img2 - tip2add
    img3 = merge(img1, img222)
    img3 = th.removeZernike(img3,zlist2adjust)
    return img3


def img2fullimg(img, fullmask, pixpos):
   # imgf = fullmask.copy()*0.
    ss = np.shape(fullmask)
    data = np.zeros(np.shape(fullmask))
    mask = np.ones(np.shape(fullmask))

    imgf = np.ma.masked_array(data,mask)
    imgsize=np.shape(img)
    imgf[pixpos[0]:pixpos[0]+imgsize[0],pixpos[1]:pixpos[1]+imgsize[1]] = img    
    return imgf

def stitch_process(coord, img_c,pixscale,zlist2fit=[1,2,3,4],zlist2adjust=[1,2,3,4]):
    fullmask=refmirr_fullmask(img_c, coord, pixscale)
    for ii in range(len(img_c[0,0])):
        img = img_c[:,:,ii]
        dm  = img2fullimg(img, fullmask , coord[ii])
        if ii == 0:
            dm0 = dm.copy()
        else:
            dm0 = coalign(dm0,dm,zlist2fit, zlist2adjust)
    return dm0

def stitch_process2(imgvec, fullmask,zlist2fit=[1,2,3,4],zlist2adjust=[1,2,3,4]):
    #img = imgvec[0]*0
    for i in range(len(imgvec)):
        img = imgvec[i]
        #dm  = img2fullimg(img, fullmask , coord[ii])
        if i == 0:
            img0 = img.copy()
        else:
            img0 = coalign(img0,img,fullmask, zlist2fit, zlist2adjust)
        print(i)
    return img0

def zernike_surf(zlist, mask, coeff=None):
    img = np.ones(np.shape(mask))
    img = np.ma.masked_array(img, mask)

    if coeff is None:
        coeff = np.ones(len(zlist))
    ccx, zmatx = zern.zernikeFitAuxmask(img, np.invert(img.mask).astype(int), zlist)

    surf = zern.zernikeSurface(img, coeff, zmatx)
    return surf

   



















 







