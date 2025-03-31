import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from arte.utils import rebin
from m4.configuration import config_folder_names as foldname
from m4.analyzers import timehistory as th
from m4.ground import geo, zernike as zern, read_data as rd, read_ottcalib_conf as roc
from m4.ground.timestamp import Timestamp
from m4.ground import read_ottcalib_conf
from m4.utils import osutils as osu
import m4.utils.parabola_footprint_registration as pr
from m4.utils.parabola_identification import ParabolaActivities

pa = ParabolaActivities()
pfr = pr.ParabolaFootprintRegistration()
OPDSERIES = foldname.OPD_SERIES_ROOT_FOLDER

def crop_frame(imgin):
    cir = geo.qpupil(-1*imgin.mask+1)
    cir = np.array(cir[0:3]).astype(int)
    img = imgin.data[cir[0]-cir[2]:cir[0]+cir[2]+1,cir[1]-cir[2]:cir[1]+cir[2]+1]
    m = imgin.mask[cir[0]-cir[2]:cir[0]+cir[2]+1,cir[1]-cir[2]:cir[1]+cir[2]+1]
    img = np.ma.masked_array(img, m)
    return img


def getMarkers(tn,flip=False, diam=24, thr=0.2):
    npix = 3.14*(diam/2)**2
    fl = osu.getFileList(tn, fold=OPDSERIES)
    nf = len(fl)
    pos = np.zeros([2,25,nf])
    for j in range(nf):
        img = th.frame(j,fl)
        if flip is True:
            print('flipping')
            img = np.fliplr(img)
        imaf = pa.rawMarkersPos(img)
        c0 = pa.filterMarkersPos(imaf, (1-thr)*npix, (1+thr)*npix)
        if j == 0:
            nmark = np.shape(c0)[1]
            pos = np.zeros([2,nmark,nf])
        pos[:,:,j]=c0
    pos = np.average(pos,2)
    return pos



def coord2ottcoord(vec1, off, flipOffset=True):
    off1 = off.copy()
    if flipOffset == True:
        off1 = np.flip(off)
        print('Offset values flipped:'+str(off1))
    vec = vec1.copy()
    for ii in range(np.shape(vec)[1]):
        vec[:,ii] = vec[:,ii] + off1
    return vec

def marker_remap(cghf,ottf):
    polycoeff = pfr.fit_trasformation_parameter(cghf,ottf)
    base_cgh = pfr._expandbase(cghf[0,:], cghf[1,:])
    cghf_tra = np.transpose(np.dot(np.transpose(base_cgh),np.transpose(polycoeff)))
    return cghf_tra

def marker_general_remap(cghf,ottf,pos2t):
    '''
        transforms the pos2t coordinates, using the cghf and ottf coordinates to create the trnasformation
    '''
    polycoeff = pfr.fit_trasformation_parameter(cghf,ottf)
    base_cgh = pfr._expandbase(pos2t[0,:], pos2t[1,:])
    cghf_tra = np.transpose(np.dot(np.transpose(base_cgh),np.transpose(polycoeff)))
    return cghf_tra


def par_remap(cgh_image, ott_image, cghf, ottf, forder=10):
    par_on_ott, mask_float, cgh_tra, difference = pfr.image_transformation(cgh_image, ott_image, cghf, ottf, forder=forder)
    cgh_tra = np.ma.masked_array(cgh_tra.data,cgh_tra == 0)
    #cgh_tra = zern.removeZernike(cgh_tra,[1,2,3,4])
    return cgh_tra 

def par_remap_only(cgh_image, cghf, ottf, forder=10):
    par_on_ott, mask_float, cgh_tra, difference = pfr.image_transformation(cgh_image, cgh_image*0, cghf, ottf, forder=forder)
    cgh_tra = np.ma.masked_array(cgh_tra.data,cgh_tra == 0)
    #cgh_tra = zern.removeZernike(cgh_tra,[1,2,3,4])
    return cgh_tra

def ott_remap(cgh_tra, ott_image):
    mask = ott_image.mask
    par_on_ott = np.ma.masked_array(cgh_tra.data,mask)
    par_on_ott = zern.removeZernike(par_on_ott,[1,2,3,4])
    ott_image = zern.removeZernike(ott_image,[1,2,3,4])
    return par_on_ott, ott_image

def marker_data_all(cgh_tn_marker,ott_tn_marker,off_cgh_marker,off_ott_marker,mark_cgh,mark_ott):
    tn0= cgh_tn_marker
    tn1 = ott_tn_marker
    fl0 = osu.getFileList(tn0, fold=OPDSERIES)
    img0=th.frame(0, fl0)
    img0= np.fliplr(img0)
    fl1 = osu.getFileList(tn1, fold=OPDSERIES)
    img1=th.frame(0, fl1)

    p0 = getMarkers(tn0,flip=True, diam=24)
    p1 = getMarkers(tn1, diam=28)

#fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,6))
#ax1.imshow(img0)
#for i in range(shape(p0)[1]):
#    ax1.text(p0[1,i],p0[0,i],i)
#ax2.imshow(img1)
#for i in range(shape(p1)[1]):
#    ax2.text(p1[1,i],p1[0,i],i)

    p0 = coord2ottcoord(p0, off_cgh_marker)
    p1 = coord2ottcoord(p1, off_ott_marker)
    pcgh = p0[:,mark_cgh]
    pott = p1[:,mark_ott]
    return pcgh,pott

def view_markers(p0, p1):
    fig, ax = plt.subplots()#figure()
    plt.plot(p0[1,:],p0[0,:],'o')
    ax.axis('equal')
    plt.plot(p1[1,:],p1[0,:],'x')
    for i in range(np.shape(p0)[1]):
        ax.text(p0[1,i],p0[0,i],str(i))        #modRB 20240518 was text(p0[1,:],p0[0,:],str(i))
        ax.text(p1[1,i]+10,p1[0,i]+10,str(i))  #text(p1[1,:]+10,p1[0,:]+10,str(i))
    plt.title('Markers position comparison')
    plt.show()
    p00 = marker_remap(p0, p1)
    dd = np.sqrt((p00[0,:]-p1[0,:])**2+(p00[1,:]-p1[1,:])**2)
    fig, ax, plt.subplots() #figure()
    plt.scatter(p1[0,:],p1[1,:],dd*50,dd)
    for i in range(np.shape(p0)[1]):
        ax.text(p1[0,i],p0[1,i],str(i))
    plt.title('Remapping error')
    plt.colorbar()

def plot_markers(p0):
    fig, ax = plt.subplots()
    plt.plot(p0[1,:],p0[0,:],'o')
    ax.axis('equal')
    plt.title('Markers position in the frame')
    for i in range(np.shape(p0)[1]):
        ax.text(p0[1,i],p0[0,i],str(i)) 

def markers_explorer(tn):
    p0=marker_data(tn, None, diam=28, flip=False)
    plot_markers(p0)
    plt.xlim(0,2048)
    plt.ylim(0,2048)
    plt.title(tn)

def marker_data(tn_marker,mark_list, diam,flip=False):
    '''
    usage:
    for ott markers: marker_data(tn,mark_list, 28,flip=False)
    for cgh markers: marker_data(tn,mark_list, 24,flip=True)
    if tn_marker is a tnvector, mark list shall be a 2D vector
    '''
    ttt = """for ott markers: marker_data(tn,mark_list, 28,flip=False)
    for cgh markers: marker_data(tn,mark_list, 24,flip=True)
    if tn_marker is a tnvector, mark list shall be a 2D vector')"""
    print(ttt)
    fl0 = osu.getFileList(tn_marker, fold=OPDSERIES)
    img0=th.frame(0, fl0)
    if flip is True:
        img0= np.fliplr(img0)
        print('flipping the frame')
    off_marker = (th.readFrameCrop(tn_marker))[2:4]
    p0 = getMarkers(tn_marker,flip, diam)
    p0 = coord2ottcoord(p0, off_marker)
    if mark_list is not None:
        p0 = p0[:,mark_list]
    return p0


def image_data(tn_img, flip=False):
    path = foldname.OPD_SERIES_ROOT_FOLDER+'/'
    fname = path + tn_img + '/average.fits'
    img = rd.read_phasemap(fname)
    if flip is True:
        img = np.fliplr(img)
    #offs = (th.readFrameCrop(tn_img))[2:4]
    conf = osu.getCameraSettings(tn_img)
    offs = [conf['x-offset'], conf['y-offset']]
    img = th.frame2ottFrame(img, offs)
    return img

def read_init_data(tnconf):
    cgh_tn_marker,cgh_tn_img,tnpar,mark_cgh_list,f0,f1,ott_tn_marker,ott_tn_img,mark_ott_list,px_ott = roc.gimmetheconf(tnconf)
    return cgh_tn_marker,ott_tn_marker,mark_cgh_list,mark_ott_list,cgh_tn_img, ott_tn_img

def init_data(tnconf):
    #(cgh_tn_marker,ott_tn_marker,mark_cgh_list,mark_ott_list,cgh_tn_img, ott_tn_img):
    '''
    '''
    cgh_tn_marker,ott_tn_marker,mark_cgh_list,mark_ott_list,cgh_tn_img, ott_tn_img = read_init_data(tnconf)
    cghf = marker_data(cgh_tn_marker,mark_cgh_list, 24,flip=True)
    #here modified
    ntn = len(ott_tn_marker)
    if ntn>1:
        print('Multi Tracknum')
        nlen=[]
        for i in mark_ott_list:
            nlen.append(len(i))
        print('N tracknum'+str(ntn))
        ottf = np.zeros([2,np.sum(nlen)])
        pos = 0
        for i in np.arange(ntn):
            ottf[:,pos:pos+nlen[i]] = marker_data(ott_tn_marker[i],mark_ott_list[i], 28,flip=False)
            pos = pos+nlen[i]

    else:
        ottf = marker_data(ott_tn_marker[0],mark_ott_list, 28,flip=False)

    '''    
    if ntn < 15:
        print('Multi Tracknum')
        nlen=[]
        for i in mark_ott_list:
            nlen.append(len(i))
        print('N tracknum'+str(ntn))
        ottf = np.zeros([2,np.sum(nlen)])
        pos = 0
        for i in np.arange(ntn):
            ottf[:,pos:pos+nlen[i]] = marker_data(ott_tn_marker[i],mark_ott_list[i], 28,flip=False)
            pos = pos+nlen[i]

    else:   
        ottf = marker_data(ott_tn_marker,mark_ott_list, 28,flip=False)
    '''
    #end of modif, keep the last line, unindented
    cgh_image = image_data(cgh_tn_img, flip=True)
    ott_image = image_data(ott_tn_img, flip=False)
    return cgh_image, ott_image, cghf, ottf

def register_par(tnconf, show=False, forder=10):
    #cgh_tn_marker,ott_tn_marker,mark_cgh_list,mark_ott_list,cgh_tn_img, ott_tn_img, show=False):
    cgh_image, ott_image, cghf, ottf = init_data(tnconf)
    if show is not False:
        view_markers(cghf, ottf)
        
    cgh_tra = par_remap(cgh_image, ott_image, cghf, ottf, forder=forder)
    cgh_tra = np.ma.masked_array(cgh_tra.data,cgh_tra == 0)
    cgh_tra = zern.removeZernike(cgh_tra,[1,2,3,4])
    ott_image = zern.removeZernike(ott_image,[1,2,3,4])
    #tn = save_registration(cgh_tra,cgh_tn_img,cgh_tn_marker,ott_tn_marker)
    tn = save_registration(cgh_tra,tnconf)

    return cgh_tra, ott_image, tn

def register_par_only(tnconf, show=False, forder=10):
    #cgh_tn_marker,ott_tn_marker,mark_cgh_list,mark_ott_list,cgh_tn_img, ott_tn_img, show=False):
    cgh_image, ott_image, cghf, ottf = init_data(tnconf)
    if show is not False:
        view_markers(cghf, ottf)

    cgh_tra = par_remap(cgh_image, ott_image, cghf, ottf, forder=forder)
    cgh_tra = np.ma.masked_array(cgh_tra.data,cgh_tra == 0)
    cgh_tra = zern.removeZernike(cgh_tra,[1,2,3,4])
    ott_image = zern.removeZernike(ott_image,[1,2,3,4])
    #tn = save_registration(cgh_tra,cgh_tn_img,cgh_tn_marker,ott_tn_marker)
    tn = save_registration(cgh_tra,tnconf)
    return cgh_tra, ott_image, tn


def save_registration(img, tnconf):#(img,cgh_tn_img,cgh_tn_marker,ott_tn_marker):
    tn = Timestamp.now()
    print(tn)
    fold = foldname.PARABOLA_REMAPPED_FOLDER+'/'+tn+'/'
    os.mkdir(fold)
    name = fold+'par_remapped.fits'
    pyfits.writeto(name, img.data)
    pyfits.append(name, img.mask.astype(int))
    copyConf(tnconf, tn)
    '''
    fobj = open(fold+'registration_info.txt','w')
    s = 'PAR CGH Tracknum: '+cgh_tn_img+'\n'
    fobj.write(s)
    s = 'PAR CGH Markers Tracknum: '+cgh_tn_marker+'\n'
    fobj.write(s)
    s = 'OTT Markers Tracknum: '+ott_tn_marker+'\n'
    fobj.write(s)
    fobj.close
    '''
    return tn
def copyConf(tnconf, tnpar):
    fromf = read_ottcalib_conf.basepath+'/'+read_ottcalib_conf.fold+tnconf+'.ini'
    tof = foldname.PARABOLA_REMAPPED_FOLDER+'/'+tnpar+'/'+tnconf+'.ini'
    shutil.copyfile(fromf, tof)

def load_registeredPar(tn,zlist = [1,2,3,4]):
    fold = foldname.PARABOLA_REMAPPED_FOLDER+'/'+tn+'/'
    name = fold+'par_remapped.fits'
    hdu = pyfits.open(name)
    img = hdu[0].data
    mask = hdu[1].data
    imgout = np.ma.masked_array(img, mask)
    imgout = zern.removeZernike(imgout, zlist)
    return imgout

def load_ott(tn, zlist=[1,2,3,4]):
    imgout=image_data(tn)
    imgout = zern.removeZernike(imgout, zlist)
    return imgout

def image_remask(img0, img1):
    img = np.ma.masked_array(img0.data, img1.mask)
    return img

def ott_calib(tnott, tnpar, zlist=[1,2,3,4,7,8], filtfreq=None, px_ott=0.00076, crpar=None):
    '''
    Inputs: tnott: tracknum OTT data in OPDSeries;
            tnpar: tracknum PAR registered, in ParabolaRemapped, produced with register_par (includes saving)
    '''
    par = load_registeredPar(tnpar)
    #modRB 20231031 to implement global fitting of Zern modes
    cir = geo.qpupil(-1*par.mask+1)
    mm = geo.draw_mask(par.data*0,cir[0],cir[1],1.44/0.00076/2,out=0)
    #mm = np.ma.masked_array(mm, mm==0)
    #par = zern.removeZernike(par, zlist) #if the par is distorted, better not fitting zernike...
    ott = load_ott(tnott)
    #ott = zern.removeZernike(ott,zlist)
    ott = zern.removeZernikeAuxMask(ott, mm, zlist)

    if filtfreq is not None:
        par = th.comp_filtered_image(par,  d=0.00076, verbose=True, disp=False, freq2filter=filtfreq)
    par = image_remask(par, ott)
    #par = zern.removeZernike(par, zlist)
    par = zern.removeZernikeAuxMask(par,mm, zlist)

    res = ott - 2*par
    #res = zern.removeZernike(res, zlist)
    res = zern.removeZernikeAuxMask(res,mm, zlist)

    view_calibration(ott, par,res, vm=50e-9,crpar=crpar)

    return res, ott, par


def mask_center(img,radius, out=0):
    print(out)
    cir = geo.qpupil(-1*img.mask+1)
    if out != 0:
        center_mask = geo.draw_mask(np.ones(img.shape),cir[0], cir[1],radius, out=1)
    else:
        center_mask = geo.draw_mask(np.zeros(img.shape),cir[0], cir[1],radius)
    
    master_mask = (-1*img.mask+1)*(-1*center_mask+1)
    master_mask = -1*master_mask+1
    img = np.ma.masked_array(img.data, master_mask)
    return img


def mask_edge(img,pix2remove):
    cir = geo.qpupil(-1*img.mask+1)
    center_mask = geo.draw_mask(np.ones(img.shape),cir[0], cir[1],cir[2]-pix2remove, out=1)
    master_mask = (-1*img.mask+1)*(-1*center_mask+1)
    master_mask = -1*master_mask+1
    img = np.ma.masked_array(img.data, master_mask)
    return img



def view_calibration(imgott, imgpar,imgres, vm=50e-9, crpar=None,nopsd=0):
    '''
        function to visualize the subtraction results.
        inputs: imgott: image of the OTT
                imgpar: registered image of the par, filtered in the case
                vm: colorbar limits
    '''
    #step1: viewing the frames
    view_subplots(imgott, imgpar, imgres, crpar=None, vm=vm, nopsd=nopsd)
    view_subplots(imgott, imgpar, imgres, crpar,vm,nopsd=nopsd)

def view_subplots(imgott, imgpar, imgres, crpar=None, vm=50e-9,nopsd=0):
    if crpar is not None:
        x = crpar[0];y=crpar[1];c=crpar[2]
        oc = imgott[x:x+c,y:y+c] ; oc = zern.removeZernike(oc,[1,2,3])
        pc = 2*imgpar[x:x+c,y:y+c]; pc = zern.removeZernike(pc,[1,2,3])
        rc = imgres[x:x+c,y:y+c];    rc = zern.removeZernike(rc,[1,2,3])
    else:
        oc = imgott.copy()
        pc = 2*imgpar
        rc = imgres.copy()

    plt.figure(figsize=(18,6))
    ax1 = plt.subplot(1,3,1)
    ax1.imshow(oc, vmin=-vm, vmax=vm)
    ax1.colorbar()
    ax1.set_title('OTT Image'+rmstitle(oc))
    
    ax2 = plt.subplot(1,3,2)
    ax2.imshow(pc, vmin=-vm, vmax=vm)
    ax2.colorbar()
    ax2.set_title('2Par Image'+rmstitle(pc))
    
    ax3 = plt.subplot(1,3,3)
    ax3.imshow(rc, vmin=-vm, vmax=vm)
    ax3.colorbar()
    ax3.set_title('Residue'+rmstitle(rc))
    
    if crpar is not None and nopsd == 0:
        norm='ortho'
        px_ott = 0.00076
        plt.clf()
        plt.figure()
        xo,yo = th.comp_psd(oc, d=px_ott, norm=norm, verbose=True)
        xp,yp = th.comp_psd(pc, d=px_ott, norm=norm, verbose=True)
        xr,yr = th.comp_psd(rc, d=px_ott, norm=norm, verbose=True)
        plt.plot(xo[1:],yo[1:]*xo[1:],'o')
        plt.plot(xp[1:],yp[1:]*xp[1:],'x')
        plt.plot(xr[1:],yr[1:]*xr[1:],'r')
        plt.legend(['OTT','2Par','Res'])
        plt.xscale('log')
        plt.yscale('log')
        plt.grid()

def rmstitle(rr):
    out=(' SfE= '+str(int(rr.std()*1e9))+'nm')
    return out

def thresh_image(img1,threshold,zlist=[1,2,3,4],inrad=0, out=0,crop=False):
    '''
    draw a central mask, as % of the radius
    '''
    #if mask_cen != 0:
    #    img = mask_center(img,mask_cen)
    img = zern.removeZernike(img1,zlist)
    idx = np.where(abs(img)>threshold)
    img.mask[idx] = True
    if inrad != 0:
        cir = geo.qpupil(-1*img.mask+1)
        img = mask_center(img,round(cir[2]*inrad))
    if out != 0:
        cir = geo.qpupil(-1*img.mask+1)
        img = mask_center(img,round(cir[2]*out), out=1)
    if crop == True:
        img = crop_frame(img) 
    return img


def quick243(img, pixs):  #pixs = [pix/m]
    img1 = crop_frame(img)
    dd = 0.03
    pp =int(pixs * dd)#.astype(int)
    ss = (np.array(img1.shape)/pp).astype(int)
#    ss = np.array((img1.shape[0]/pp)).astype(int)
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

def quick283(img, pixs):  #pixs = [pix/m]
    img1 = crop_frame(img)
    dd = 0.08
    pp =int(pixs * dd)#.astype(int)
    ss = (np.array(img1.shape)/pp).astype(int)
#    ss = np.array((img1.shape[0]/pp)).astype(int)
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

def adjust_marker(cghf,ottf, mid,ran):
    cf=marker_remap(cghf,ottf)
    dd0 = np.sqrt((cf[0,:]-ottf[0,:])**2+(cf[1,:]-ottf[1,:])**2)
    x=np.linspace(-ran,ran,2*ran)
    y=np.linspace(-ran,ran,2*ran)
    w = np.zeros((2*ran,2*ran))
    for ii in np.arange(len(x)):
        for jj in np.arange(len(y)):
            tmp = cghf.copy()
            tmp[:,mid]=tmp[:,mid]+[x[ii],y[jj]]
            cf=marker_remap(tmp,ottf)
            dd = np.sqrt((cf[0,:]-ottf[0,:])**2+(cf[1,:]-ottf[1,:])**2)
            w[ii,jj] = dd.std()
    a= np.array(np.where(w==np.min(w))).flatten()
    plt.imshow(w)
    plt.colorbar()
    plt.title('Marker scatter')
    print(a)
    ss = w.shape
    offs = np.array([x[a[0]],y[a[1]]])
    print(offs)
    tmp = cghf.copy()
    tmp[:,mid]=tmp[:,mid]+offs
    print('Initial pos error sigma:');print(dd0.std())
    print('Minimum pos error sigma:');print(np.min(w))
    return tmp

def compSlope(img,px, rfact):
    ss=np.array(np.shape(img))
    sli = img-np.roll(img,(1,1),axis=(0,1))
    slid = geo.congrid2D(sli,(ss/rfact).astype(int))/(px*rfact)
    slim = geo.congrid2D(-1*sli.mask+1,(ss/rfact).astype(int))
    sli = np.ma.masked_array(slid,(-1*slim+1))
    return sli.std()

def compSlopXY(img,px, rfact):
    ss=np.array(np.shape(img))
    slm  = -1*img.mask+1
    slm = slm+np.roll(slm,(1,1),axis=(0,1))
    slix = np.ma.masked_array((img-np.roll(img,(1,0),axis=(0,1)))/px,slm < 2)
    sliy = np.ma.masked_array((img-np.roll(img,(0,1),axis=(0,1)))/px,slm < 2)
    sli  = np.ma.sqrt(slix**2+sliy**2)
    slid = geo.congrid2D(sli,(ss/rfact).astype(int))
    slim = geo.congrid2D(-1*sli.mask+1,(ss/rfact).astype(int))
    sli = np.ma.masked_array(slid,(-1*slim+1))
    return sli.std()


def compSlopXY2(img,px, rfact, thr=None):
    ss=np.array(np.shape(img))
    imgr = rebin.rebin(img, (ss/rfact).astype(int))
    slm = -1*imgr.mask+1
    #imgr = geo.congrid2D(img,(ss/rfact).astype(int))
    #imgm = geo.congrid2D(-1*img.mask+1,(ss/rfact).astype(int))
    #imgr = np.ma.masked_array(imgr,(-1*imgm+1))
    #slm  = imgm
    slm = slm+np.roll(slm,(1,1),axis=(0,1))
    #slix = np.ma.diff(imgr,1,axis=0)/(px*rfact)
    #sliy = np.ma.diff(imgr,1,axis=1)/(px*rfact)
    slix = np.ma.masked_array((imgr-np.roll(imgr,(1,0),axis=(0,1)))/(px*rfact),slm < 2)
    sliy = np.ma.masked_array((imgr-np.roll(imgr,(0,1),axis=(0,1)))/(px*rfact),slm < 2)
    sli  = np.ma.sqrt(slix**2+sliy**2)
    if thr is not None:
        sli = np.ma.masked_array(sli, sli.data > thr)
    #sli = np.ma.masked_array(sli, slm < 2)
    #sli = np.ma.masked_array(slid,(-1*slim+1))
    return sli
