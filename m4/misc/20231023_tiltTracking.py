import numpy as np
from matplotlib import pyplot as plt
from m4.misc import par_meas as pp
from m4.ground import read_data as rr
from m4.mini_OTT import timehistory as th
from m4.ground import zernike as zern
from m4.mini_OTT.measurements import Measurements
from m4 import main
from astropy.io import fits as pyfits
from m4.configuration import start
from m4.devices.i4d import I4D
from m4.configuration.ott_parameters import Interferometer
from m4 import noise
import time
from m4.misc import image_registration_lib as imgreg
from m4.ground import geo

conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
ott, interf, dm = start.create_ott(conf)
meas = Measurements(ott,interf)

tnc = '20230720_192515'
zern2corrf = np.array([2]);dofidf = np.array([0])# focus, parpist
zern2corrf = np.array([0,1,3,4]); dofidf = np.array([1,2,3,4])# TipTilt, Coma,ParTip, ParTilt, RmTip, RmTilt
zern2corrf = np.array([0,1]) ;dofidf = np.array([3,4])#TipTilt, RmTip, RmTilt


#acquisition
rmstep = 40
deltapos=540
nstep = 20
zern2corrf = np.array([0,1]) ;dofidf = np.array([3,4])#TipTilt, RmTip, RmTilt
rm0= ott.referenceMirrorSlider.getPosition()
tnlist =[]
for i in np.arange(nstep):
    ott.referenceMirrorSlider.setPosition(rm0+rmstep*i)
    time.sleep(10)
    tna = main.align_PARAndRM(ott, interf, tnc, zern2corrf, dofidf, n_frames=16, doit=True)
    tnlist.append(tna)

#the following list if from 130 to 890, step 40. aligned coma and power at 130
tnlist = ['20231024_155008', '20231024_155130', '20231024_155250', '20231024_155409', '20231024_155532', '20231024_155654', '20231024_155824', '20231024_155947', '20231024_160111', '20231024_160234', '20231024_160356', '20231024_160518', '20231024_160641', '20231024_160806', '20231024_160927', '20231024_161050', '20231024_161214', '20231024_161337', '20231024_161501', '20231024_161625']

#second run, anitial alignment at rm=600 with 0coma and power
tnlist = ['20231025_131146', '20231025_131306', '20231025_131427', '20231025_131546', '20231025_131706', '20231025_131825', '20231025_131946', '20231025_132107', '20231025_132228', '20231025_132348', '20231025_132509']

#third run, with -22 nm coma at center
tnlist = ['20231025_134243', '20231025_134400', '20231025_134520', '20231025_134639', '20231025_134759', '20231025_134921', '20231025_135041', '20231025_135202', '20231025_135323', '20231025_135444', '20231025_135605']

#4th run, starting from center, global alignment implemented, 12 nm global coma at center, par not removed.
tnlist = ['20231027_113212', '20231027_113334', '20231027_113459', '20231027_113624', '20231027_113749', '20231027_113914', '20231027_114041', '20231027_114208', '20231027_114334', '20231027_114458', '20231027_114623']

#5th run, starting from center, global alignment implemented, 13 nm global coma at center, PAR removed
tnlist = ['20231027_123944', '20231027_124106', '20231027_124231', '20231027_124355', '20231027_124519', '20231027_124643', '20231027_124807', '20231027_124932', '20231027_125056', '20231027_125220', '20231027_125345']

#analysis
#testing with auxilimary mask def zernikeFitAuxmask(img, auxmask, zernike_index_vector):

tnpar  = '20231016_124531'
par = imgreg.load_registeredPar(tnpar) 
cir = geo.qpupil(-1*par.mask+1)
mm = geo.draw_mask(par.data*0,cir[0],cir[1],1.44/0.00076/2,out=0)
tt = []
zz = []
zg = []
zr = []
zrg= []
b = '/mnt/m4storage/Data/M4Data/OPTData/Alignment/'
for i in tnlist:
    h=(pyfits.open(b+i+'/PositionAndDeltaCommand.fits'))[0].data
    tt.append(h[1,3:5])
    q=th.read_phasemap(b+i+'/StartImage.fits')
    #q=th.read_phasemap(b+i+'/FinalImage.fits')

    q = th.frame2ottFrame(q,[580,20])
    cc, m=th.zernike.zernikeFit(q,[1,2,3,4,5,6,7,8])
    cc1, m = zern.zernikeFitAuxmask(q,mm,[1,2,3,4,5,6,7,8])
    res = q-2*par
    cr, m=th.zernike.zernikeFit(res,[1,2,3,4,5,6,7,8])
    cr1, m = zern.zernikeFitAuxmask(res,mm,[1,2,3,4,5,6,7,8])
    zz.append(cc)
    zg.append(cc1)
    zr.append(cr)
    zrg.append(cr1)
    

tt=np.array(tt)
zz = np.array(zz)
zg = np.array(zg)
zr = np.array(zr)
zrg = np.array(zrg)


figure();grid('on')
plot(zz[:,3],'-o')
plot(zz[:,4],'-o')
plot(zz[:,5],'-o')
plot(zz[:,6],'-o')
plot(zz[:,7],'-o')
legend(['Z4','Z5','Z6','Z7','Z8']);title('Local ZModes, from ParCenter upward')

figure();grid('on')
plot(zg[:,3],'-o')
plot(zg[:,4],'-o')
plot(zg[:,5],'-o')
plot(zg[:,6],'-o')
plot(zg[:,7],'-o')
legend(['Z4','Z5','Z6','Z7','Z8']);title('Global ZModes, from ParCenter upward')

figure();grid('on')
plot(zr[:,3],'-o')
plot(zr[:,4],'-o')
plot(zr[:,5],'-o')
plot(zr[:,6],'-o')
plot(zr[:,7],'-o')
legend(['Z4','Z5','Z6','Z7','Z8']);title('Local ZModes on Res, from ParCenter upward')

figure();grid('on')
plot(zrg[:,3],'-o')
plot(zrg[:,4],'-o')
plot(zrg[:,5],'-o')
plot(zrg[:,6],'-o')
plot(zrg[:,7],'-o')
legend(['Z4','Z5','Z6','Z7','Z8']);title('Global ZModes on Res, from ParCenter upward')

