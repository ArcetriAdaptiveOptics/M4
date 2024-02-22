'''
this is the master script for the OTT map calibration and REQ verification.
data are not collected with this script. This script is just for data analysis.

the workflow for the OTT calibration and  verification is a follows:
-preparation of a calibration package:
-- in image_registration lib
    
usage of the calibration
-- in ott_calibration_lib

REQ verification
 in the current script

how to add markers in 20230901_par_registration_test.py
'''
import numpy as np
from m4.utils import image_registration_lib as imgreg
from m4.mini_OTT import timehistory as th
from m4.ground import read_ottcalib_conf
from m4.ground import zernike as zern
from m4.analyzers import requirement_analyzer as ra
from m4.ground import geo
from matplotlib.pyplot import *
from importlib import reload
from astropy.io import fits as pyfits
# Initilization data

px_ott = 0.00076
f0 = 0
f1 = 40

tnconf = '20231011_100000' #old L3. Configuration to remap the PAR
tnconf = '20231011_140000' #new L3. Configuration to remap the PAR
tnpar  = '20231011_150952' #associated ParRemapped
tnconf = '20231013_230000' # 3 tracknum for ott markers
tnpar  = '20231016_121641' #order 6 '20231014_005211' #associated ParRemapped
tnpar  = '20231016_124531' #order 10

tnott  = '20231011_164520' #RM270
tnott  = '20231012_224214' #RM450
tnott  = '20231006_123231' #RM600
tnott  = '20231007_114524' #RM600
tnott  = '20231012_182010' #RM750
tnott  = '20231010_210434' #RM930

#OTT meas after alignment on local coordinates
tnott  = '20231011_164520' #RM270
tnott  = '20231012_224214' #RM450
tnott  = '20231006_123231' #RM600
tnott  = '20231007_114524' #RM600
tnott  = '20231012_182010' #RM750
tnott  = '20231010_210434' #RM930
#----
#OTT meas after alignment on global coordinates. tnpar used=20231016_124531
tnott = '20231029_150436' #RM200
tnott = '20231029_105812' #RM400
tnott = '20231029_175926' #RM600
tnott = '20231028_223658' #RM800
tnott = '20231028_204448' #RM1000
z2corr = [1,2,3,4]


#### HOW TO SECTION
#restoring the configuration (All TNs and markers list)
cgh_tn_marker,cgh_tn_img,tnpar,mark_cgh_list,f0,f1,ott_tn_marker,ott_tn_img,mark_ott_list,px_ott = read_ottcalib_conf.gimmetheconf(tnconf)

#How to create a registered image of the PAR
#processing the data
par_remapped, ott_image, tnpar = imgreg.register_par(tnconf)
#or reloading the data
par_remapped = imgreg.load_registeredPar(tnpar)
ott = imgreg.load_ott(tnott)

#how to filter the PAR and compute the OTT calibrated SFMap
par_filtered = th.comp_filtered_image(par_remapped,  d=px_ott, verbose=True, disp=False, freq2filter=(f0,f1))
par = imgreg.image_remask(par_filtered, ott)
#res = ott-2*par_ott
res = ott-2*par

#OTT calibrated SFmap in a single shot (requires the TNPAR processed)
#---> the data output here are the calibration results, to be analyzed
res, ott, par= imgreg.ott_calib(tnott, tnpar, zlist=[1,2,3,4],filtfreq=[0,40], crpar=None)
imgreg.view_calibration(ott, par, vm=50e-9,crpar=[700,850,150])



#data analysis
#analysis of the astigmatism vs RM pos (truss is fixed)
tnlist= ['20231011_164520','20231012_224214','20231006_123231', '20231007_114524','20231012_182010', '20231010_210434']
rmpos = [270,450,600,600,750,930]
ax=[]
ay = []
for i in tnlist:
    ott = imgreg.load_ott(i, zlist=[1])
    cc, m=th.zernike.zernikeFit(ott,[1,2,3,4,5,6])
    ax.append(cc[4])
    ay.append(cc[5])
figure()
plot(rmpos, ax*1e9,'o');xlabel('RM position [mm]');ylabel('Astigmatism Coeff [nm]')
plot(rmpos, ay*1e9,'o')
legend(['AstX','AstY'])


###
#analysis of the astigmatism vs Truss pos (RM is fixed)
tnlist= []
rmpos = []
ax=[]
ay = []
for i in tnlist:
    ott = imgreg.load_ott(i, zlist=[1])
    cc, m=th.zernike.zernikeFit(ott,[1,2,3,4,5,6])
    ax.append(cc[4])
    ay.append(cc[5])
figure()
plot(rmpos, ax*1e9,'o');xlabel('RM position [mm]');ylabel('Astigmatism Coeff [nm]')
plot(rmpos, ay*1e9,'o')
legend(['AstX','AstY'])


#Zernike analysis
zlist=np.arange(10)+1
zott, m = zern.zernikeFit(ott_image,zlist)
zres, m = zern.zernikeFit(res,zlist)
zparott,m = zern.zernikeFit(par_ott,zlist)


#comments
'''
1: par_on_ott ha coma, possibly due to remasking as OTT.
what is the zero coma value?
- suggestion, use auxmask == fitting the zernike on the 1.5m diameter circle
- or remove everythin in the final image
- or align the OTT to the coma value measured on the par_on_ott ones

'''

#verification of remapped markers in PAR
#to check for deformations due to partial mapping
tnconf = '20231011_140000' #new L3. Configuration to remap the PAR
cgh_image, ott_image, cghf, ottf = imgreg.init_data(tnconf)
cgh_tn_marker,cgh_tn_img,tnpar,mark_cgh_list,f0,f1,ott_tn_marker,ott_tn_img,mark_ott_list,px_ott = read_ottcalib_conf.gimmetheconf(tnconf)
parmark = imgreg.marker_data(cgh_tn_marker,None,24, flip=True)
p00 = parmark
p0 = imgreg.marker_general_remap(cghf,ottf,parmark)
plot(p0[0,:],p0[1,:],'x')
for i in range(25):
    text(p0[0,i],p0[1,i],i)
from m4.utils.parabola_identification import ParabolaActivities
pa = ParabolaActivities()
circ0 = [8,11,12]
circ1 = [6,7,13,14,17,18]
circ2 = [3,4,5,15,16,20,21]
circ3 = [0,1,2,9,10,19,22,23,24]
radii = np.array([0.1, 0.2505, 0.452, 0.6911])
c0, axs0, r0 = pa._fitEllipse(p0[0,circ0],p0[1,circ0]); c0 = np.real(c0)
c1, axs1, r1 = pa._fitEllipse(p0[0,circ1],p0[1,circ1]); c1 = np.real(c1)
c2, axs2, r2 = pa._fitEllipse(p0[0,circ2],p0[1,circ2]); c2 = np.real(c2)
c3, axs3, r3 = pa._fitEllipse(p0[0,circ3],p0[1,circ3]); c3 = np.real(c3)
pradii = np.array([r0, r1, r2, r3])
ps = pradii/radii

c0, axs0, r0 = pa._fitEllipse(p00[0,circ0],p00[1,circ0]); c0 = np.real(c0)
c1, axs1, r1 = pa._fitEllipse(p00[0,circ1],p00[1,circ1]); c1 = np.real(c1)
c2, axs2, r2 = pa._fitEllipse(p00[0,circ2],p00[1,circ2]); c2 = np.real(c2)
c3, axs3, r3 = pa._fitEllipse(p00[0,circ3],p00[1,circ3]); c3 = np.real(c3)
pradii0 = np.array([r0, r1, r2, r3])
ps0 = pradii0/radii

#***************
#again on marker check: 
tnconf0 = '20231011_140000' #new L3. Configuration to remap the PAR
tnconf1 = '20231013_230000'
cgh_tn_marker = '20230428_093328'
parmark = imgreg.marker_data(cgh_tn_marker,None,24, flip=True)

pp, oo, pf0, of0= imgreg.init_data(tnconf0)
p0 = imgreg.marker_general_remap(pf0,of0,parmark)
pp, oo, pf1, of1= imgreg.init_data(tnconf1)
p1 = imgreg.marker_general_remap(pf1,of1,parmark)
#d0 = of1-p0
#d1 = of1-pf1

plot(of1[0,:],of1[1,:],'o')
plot(of0[0,:],of0[1,:],'d')
plot(p0[0,:],p0[1,:],'+')
plot(p1[0,:],p1[1,:],'*')
axis('scaled')
legend(['OTT all mark','OTT mark center','PARmarkRemap, centerRef','ParmarkRemap, 3posRef'])

#showing only the used markers
plot(of1[0,:],of1[1,:],'o')
plot(p1[0,:],p1[1,:],'*')
axis('scaled')
legend(['OTT all markers','PAR markers Remapped'])

#load OTT

#load PAR

#*******************************************
#*** MASTER CALIBRATION AND REQ VERIFICATION
px_ott = 0.00076
res, ott, par= imgreg.ott_calib(tnott, tnpar, zlist=[1,2,3,4],filtfreq=[0,40], crpar=None)

#slope error < 0.35 as
rfact = 4
slo = imgreg.compSlope(res, px_ott, rfact);print('Slope error, pxscale: '+px_ott*rfact);print(slo*206265)

#HF error < 9 nm, with mild manipulation
res1 = imgreg.thresh_image(res, 50e-9,zlist,inrad=0.2)
qq = imgreg.quick243(res1*2, 1/px_ott);figure();imshow(qq);colorbar();title('WFE at interact scale')
qv = (np.sort(qq.flatten()));figure();plot(qv)
qs = qv[int(np.sum(np.invert(qq.mask))*0.95)]
print('HF WF error:')
print(qs)
ps = 1/0.00076
step =10000
n_patches = 3
slope = ra.test242(res1, 1/px_ott, display=False)

#roc = ra.test283(res, 1/px_ott, step, display=False) questo è da correggere con v2

diff_piston = ra.diffPiston(res*2);print('Diff piston error: '+diff_piston)
rms31 = ra.test243(res*2, 0.015, 1/px_ott, step, n_patches);print('RMS HF '+rms31)
rms500 = ra.test243(res*2, 0.25, 1/px_ott, step, n_patches);print('RMS 500 mm :'+ rms500)

#some images
rcParams['image.cmap'] = 'hot'

#images with the PAtches and the par subtraction
tnott = '20231028_204448'
tnpar = '20231016_124531'
res, ott, par= imgreg.ott_calib(tnott, tnpar, zlist=[1,2,3,4,7,8],filtfreq=None, crpar=None)
imgreg.view_calibration(ott, par, res,vm=50e-9,crpar=[1600,800,300])


base = '/mnt/backup/Archeopterix_20230517/Data/Results/Req/'
tnlist  = ['20231011_164520' , '20231012_224214', '20231006_123231','20231012_182010', '20231010_210434']
pos = [270,450,600,750,930]
rr =[]
ss = []
dp=[]
rfact = 4
doshow = 0
for i in range(len(tnlist)):
    res, ott, par= imgreg.ott_calib(tnlist[i], tnpar, zlist=[1,2,3,4,7,8],filtfreq=[0,20], crpar=None)
    dr=str(int(ott.std()*2*1e9))
    if doshow == 1:
        figure();clf();imshow(2*ott);title(tnlist[i]+' OTT image, WFE= '+dr+'nm');colorbar()
        savefig(base+tnlist[i]+'-OTT.png')
    dr=str(int(res.std()*2*1e9))
    if doshow ==1:
        figure();clf();imshow(2*res);title(tnlist[i]+' OTT calibr. image, WFE= '+dr+'nm');colorbar()
        savefig(base+tnlist[i]+'-OTTCalib.png')
    g=th.removeZernike(2*res,[1,2,3,4,5,6,7,8])
    rr.append(g.std())
    slo = imgreg.compSlope(res, px_ott, rfact)
    ss.append(slo*206265)
    dp.append( ra.diffPiston(res*2))

dp=np.array(dp)
print('diff pist');print(dp.mean);print(dp.std())
rr=np.array(rr)
ss=np.array(ss)
close('all')
m = []
for i in tnlist:
    q=imgreg.load_ott(i)
    m.append(-1*q.mask+1)
imshow(sum(m,0))

rr=[]
for i in tnlist:
    res, ott, par= imgreg.ott_calib(i, tnpar, zlist=[1,2,3,4,7,8],filtfreq=[0,18], crpar=None);roc = ra.test283(res,1/px_ott, step);print(roc)
    rr.append(roc)



#images for the report
#ratio rms to PtV for Zernike modes 2-11: PtV/RMS
zamp = np.array([20,20,20,20,20,6])
zampptv =np.array([4,2*np.sqrt(3),np.sqrt(6),np.sqrt(8),2*np.sqrt(5),3.35])
#
rcParams['image.cmap'] = 'hot'
base = '/mnt/m4storage/Data/Results/Req/'
base = '/mnt/backup/Archeopterix_20230517/Data/Results/Req/'
tnlist  =[ '20231029_150436', '20231029_105812','20231029_175926','20231028_223658', '20231028_204448']
zoomy=[200,600,900,1400,1600];zoomx=[800,800,800,900,800]
pos = [200,400,600,800,1000]
ss = []
dp=[]
rc = []
allcurv = []
rfact = 4
doshow = 0
step = 10000
vm=100e-9
cpatch= 60 #(15 patches on the RM diam)
outp= 2
cmaps = []
rfold = '/mnt/backup/Archeopterix_20230517/Data/M4Data/Result/Data4ESO/'

for i in range(len(tnlist)):
    res, ott, par= imgreg.ott_calib(tnlist[i], tnpar, zlist=[1,2,3,4,7,8],filtfreq=[0,40], crpar=None)
    if doshow == 1:
        savefig(base+tnlist[i]+'-comparison.png')
        imgreg.view_calibration(ott, par, res,vm=50e-9,crpar=[zoomy[i],zoomx[i],300], nopsd=1)
        savefig(base+tnlist[i]+'-zoom.png')
    fname = tnlist[i]+'OTT-Cal.fits';pyfits.writeto(rfold+fname, res.data);pyfits.append(rfold+fname, res.mask.astype(np.uint8))
    dr=str(int(ott.std()*2*1e9))
    if doshow == 1:
        figure();clf();imshow(2*ott,vmin=-vm, vmax=vm);title(tnlist[i]+' OTT, WFE= '+dr+'nm');colorbar()
        savefig(base+tnlist[i]+'-OTT.png')
    dr=str(int(res.std()*2*1e9))
    if doshow ==1:
        figure();clf();imshow(2*res,vmin=-vm, vmax=vm);title(tnlist[i]+' OTT calibr., WFE= '+dr+'nm');colorbar()
        savefig(base+tnlist[i]+'-OTTCalib.png')
    slo = imgreg.compSlopXY2(res,px_ott,2, 10e-6)#slo = imgreg.compSlope(res, px_ott, rfact)
    ss.append(slo.std()*206265)
    if doshow == 1:
        imshow(slo*206265);colorbar();title(tnlist[i]+' Surf Slope '+f"{slo.std()*206265:.{2}f}")
        savefig(base+tnlist[i]+'-slope.png')
    dp.append( ra.diffPiston(res*2))
    w=ra.patches_analysis_map(res, 0.04, 1/px_ott, cpatch)
    cmaps.append(w)
    w1=imgreg.mask_edge(w, outp)
    if doshow == 1:
        figure()
        imshow(w1, vmax=80);colorbar();title('Local Curvature [km]')
        savefig(base+tnlist[i]+'-Curvature.png')
    wl = w1[w1.mask == False]
    wl = wl.data
    wl = np.sort(wl.flatten())
    roc = wl[len(wl)*int(0.05)]
    rc.append(roc)
    allcurv.append(wl)

dp=np.array(dp)
rr=np.array(rr)
ss=np.array(ss)
rc = np.array(rc)

roc = np.sort(np.concatenate(allcurv))
plot(roc);title('Local curvature distribution');xlabel('Patch id');ylabel('Curvature [km]')
rx = roc[int(len(roc)*0.05)]
plot(int(len(roc)*0.05),rx,'or')
print(rx)

#hacking the curvature data
fig = figure(figsize=(20,4))
outp=4
clf()
cmask = cmaps[0].data*0
cmask[28:33,27:33 ]=1;cmask[39:45,22:28]=1;
cc = []
for i in range(0,5):
    w=cmaps[i]
    w1=imgreg.mask_edge(w, outp)
    #ow1 = w1*cmask
    mask = w1.mask
    mask[cmask==1]=True
    wl = np.ma.masked_array(w1,mask)
    subplot(1,5,i+1);imshow(wl, vmax=80)
    wl = wl[wl.mask == False]
    wl = wl.data
    wl = np.sort(wl.flatten())
    cc.append(wl)
roc = np.sort(np.concatenate(cc))
plot(roc);title('Local curvature distribution');xlabel('Patch id');ylabel('Curvature [km]')
plot(int(len(roc)*0.05),rx,'or')
rx = roc[int(len(roc)*0.05)]
print(rx)

#filtering below 20 km
clf()
cc=[]
for i in range(0,5):
    w=cmaps[i]
    w1=imgreg.mask_edge(w, outp)
    #w1 = w1*cmask
    #wl = w1[w1.mask == False]
    wl = np.ma.masked_array(w1.data*(w1.mask==0),mask=((w1.data)*(w1.mask==0)<30))
    subplot(1,5,i+1);imshow(wl, vmax=80)
    wl = wl[wl.mask ==  False]
    wl = np.sort(wl.flatten())
    cc.append(wl)
roc = np.sort(np.concatenate(cc))
rx = roc[int(len(roc)*0.05)]
print(rx)



#curvature test on flat mirror
q=th.read_phasemap('/home/labot/Desktop/FlatMirror.fits')


print('diff pist');print(dp.mean());print(dp.std())
print('WFE 8ZRemoved');print(rr.mean());print(rr.std())
print('Slope error');print(ss.mean());print(rr.std())
print('Curv error');print(rc.mean());print(rc)

def congrid(img, rfact):
    ss = np.array(np.shape(img))
    slid = geo.congrid2D(img,(ss/rfact).astype(int))
    slim = geo.congrid2D(-1*img.mask+1,(ss/rfact).astype(int))
    sli = np.ma.masked_array(slid,(-1*slim+1))
    return sli


#test about residual coma
from m4.ground import geo
N=2048
im = np.ones([N,N])
mask=geo.draw_mask(im,N/2,N/2,1895/2,out=1)
im=np.ma.masked_array(im,mask)
cc, mat = zern.zernikeFit(im,[1,2,3,4,5,6,7,8,9,10])
c6=10e-9
c7=c6
im[im.mask == False]=c6*mat[:,6]+c7*mat[:,7]



'''
rfold = '/mnt/backup/Archeopterix_20230517/Data/M4Data/Result/Data4ESO/'
fname = rfold+cgh_tn_marker+'cgh-markers-coord.fits';pyfits.writeto(fname, cghf);
fname = rfold+ott_tn_marker+'ott-markers-coord.fits';pyfits.writeto(fname, ottf);
fname = rfold+cgh_tn_img+'-PARcgh-image.fits';pyfits.writeto(fname, cgh_image.data);pyfits.append(fname, cgh_image.mask.astype(np.uint8))
fname = rfold+ott_tn_img+'-OTT-image.fits';pyfits.writeto(fname, ott_image.data);pyfits.append(fname,ott_image.mask.astype(np.uint8))
'''
