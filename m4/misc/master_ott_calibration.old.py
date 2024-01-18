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
from m4.misc import image_registration_lib as imgreg
from m4.mini_OTT import timehistory as th
from m4.misc import read_ottcalib_conf
from m4.ground import zernike as zern
from m4.analyzers import requirement_analyzer as ra
from importlib import reload
# Initilization data

px_ott = 0.00076
f0 = 0
f1 = 40

tnconf = '20231011_100000' #old L3. Configuration to remap the PAR
tnconf = '20231011_140000' #new L3. Configuration to remap the PAR
tnpar  = '20231011_150952' #associated ParRemapped
tnconf = '20231013_230000' # 3 tracknum for ott markers
tnpar  = '20231016_121641' #order 6 '20231014_005211' #associated ParRemapped
tnpar  = '20231016_124531'
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

plot(p0[0,:],p0[1,:],'x')
plot(p1[0,:],p1[1,:],'o')
plot(of0[0,:],of0[1,:],'*')
plot(of1[0,:],of1[1,:],'+')

legend(['Ref at center','Ref in 3 pos','OTT mark','OTT mark 3 pos'])
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
qv = (sort(qq.flatten()));figure();plot(qv)
qs = qv[int(np.sum(np.invert(qq.mask))*0.95)]
print('HF WF error:')
print(qs)
ps = 1/0.00076
step =10000
n_patches = 3
slope = ra.test242(res1, 1/px_ott, display=False)

roc = ra.test283(res, 1/px_ott, step, display=False)# patch analysis already uses the curv_fit_v2 code by Luca to fix the curv computation questo è da correggere con v2

diff_piston = ra.diffPiston(res*2);print('Diff piston error: '+diff_piston)
rms31 = ra.test243(res*2, 0.015, 1/px_ott, step, n_patches);print('RMS HF '+rms31)
rms500 = ra.test243(res*2, 0.25, 1/px_ott, step, n_patches);print('RMS 500 mm :'+ rms500)

#some images
tn  = ['20231011_164520' , '20231012_224214', '20231006_123231','20231012_182010', '20231010_210434']
pos = [270,450,600,750,930]
m = []
for i in tn:
    q=imgreg.load_ott(i)
    m.append(-1*q.mask+1)
imshow(sum(m,0))

rr=[]
for i in tn:    
    res, ott, par= imgreg.ott_calib(i, tnpar, zlist=[1,2,3,4,7,8],filtfreq=[0,18], crpar=None);roc = ra.test283(res,1/px_ott, step);print(roc)
    rr.append(roc)

#some more images, new meas with global alignment
tn=[ '20231029_150436', '20231029_105812','','20231028_223658', '20231028_204448']

pos = [200, 400, 600, 800, 1000]

#images for the report
rcParams['image.cmap'] = 'hot'
base = '/mnt/m4storage/Data/Results/Req/'
tnlist  =[ '20231029_150436', '20231029_105812','20231029_175926','20231028_223658', '20231028_204448']
pos = [200,400,600,800,1000]
rr =[]
ss = []
dp=[]
rc = []
rfact = 4
doshow = 1
for i in range(len(tnlist)):
    res, ott, par= imgreg.ott_calib(tnlist[i], tnpar, zlist=[1,2,3,4],filtfreq=[0,40], crpar=None)
    dr=str(int(ott.std()*2*1e9))
    if doshow == 1:
        figure();clf();imshow(2*ott);title(tnlist[i]+' OTT, WFE= '+dr+'nm');colorbar()
        savefig(base+tnlist[i]+'-OTT.png')
    dr=str(int(res.std()*2*1e9))
    if doshow ==1:
        figure();clf();imshow(2*res);title(tnlist[i]+' OTT calibr., WFE= '+dr+'nm');colorbar()
        savefig(base+tnlist[i]+'-OTTCalib.png')
    g=th.removeZernike(2*res,[1,2,3,4,5,6,7,8])
    rr.append(g.std())
    slo = imgreg.compSlope(res, px_ott, rfact)
    ss.append(slo*206265)
    dp.append( ra.diffPiston(res*2))
    roc = ra.test283(res,1/px_ott, step)
    rc.append(roc)

dp=np.array(dp)
print('diff pist');print(dp.mean());print(dp.std())
rr=np.array(rr)
print('WFE 8ZRemoved');print(rr.mean());print(rr.std())
ss=np.array(ss)
print('Slope error');print(ss.mean());print(rr.std())
rc = np.array(rc)

