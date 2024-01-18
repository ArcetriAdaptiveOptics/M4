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
conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
ott, interf, dm = start.create_ott(conf)
meas = Measurements(ott,interf)
phcamfocus = I4D(Interferometer.i4d_IP, Interferometer.i4d_port)
#  4D Configurations

aconf = 'D:/config/20231002_ott_RM975-900x1100_600-880_AlignMask_aveRM.4Dini'
mconf = 'D:/config/20231002_ott_RM975-900x1100_600-880_noMask_aveRM.4Dini'

tnc = '20230720_192515'
zern2corrf = np.array([0,1]) ;dofidf = np.array([3,4])#TipTilt, RmTip, RmTilt

###*** moving the RM carriage
rmstep = -10
startpos = 1000
endpos = 300
nstep = int(abs((endpos-startpos)/rmstep))
rms0=ott.referenceMirrorSlider.getPosition()
tnlist = []

for i in range(nstep):
    print('Step id:'+str(i))
    print('Moving the RM to: '+str(rms0+i*rmstep))
    time.sleep(3)
    ott.referenceMirrorSlider.setPosition(rms0+rmstep*i)
    time.sleep(10)
    #phcamfocus.loadConfiguration(aconf)
    tna = main.align_PARAndRM(ott, interf, tnc, zern2corrf, dofidf, n_frames=4, doit=True)
    #phcamfocus.loadConfiguration(mconf)
    tn=meas.opticalMonitoring(100,3)
    tnlist.append(tn)


###moving the Truss carriage
ttstep = -10
startpos = 880
endpos = 180
nstep = int(abs((endpos-startpos)/ttstep))
tt0=ott.parabolaSlider.getPosition()
tnlist = []

for i in range(nstep):
    print('Step id:'+str(i))
    print('Moving the PAR to: '+str(tt0+i*ttstep))
    time.sleep(3)
    ott.parabolaSlider.setPosition(tt0+ttstep*i)
    time.sleep(10)
    #phcamfocus.loadConfiguration(aconf)
    tna = main.align_PARAndRM(ott, interf, tnc, zern2corrf, dofidf, n_frames=4, doit=True)
    #phcamfocus.loadConfiguration(mconf)
    tn=meas.opticalMonitoring(100,3)
    tnlist.append(tn)


tnlist=['20231017_090413','20231017_090639','20231017_090805','20231017_090955','20231017_091054','20231017_091339','20231017_091454','20231017_091556','20231017_091651','20231017_091750']
tnlist = ['20231017_0904130','20231017_091651','20231017_091556','20231017_091454','20231017_090639','20231017_091339','20231017_090805','20231017_090955']
#tnlist = [
ttpos = [880,480,280,180,100,380,580,680,780,880]
tn1= '
#data analysis
tnlist =['20231002_130819','20231002_131350','20231002_131931','20231002_132511','20231002_133043','20231002_133621','20231002_134203','20231002_134743','20231002_135324','20231002_135859','20231002_140444',]
tn0='20231002_130819'
tn1 = '20231002_191625'
tnlist = th.tnscan(tn0, tn1)

ntn = len(tnlist)
qc = []
for i in arange(ntn):
    fl = th.fileList(tnlist[i])
    q = th.averageFrames(0,99,fl)
    q = th.removeZernike(q,[1,2,3,4,5,6,7,8,9,10])
    qc.append(q) 

x =400;y = 300;c=200; vm=50e-9
qq=np.ma.masked_array(qc)
z=np.ma.mean(qq,axis=0)
z1=z[x:x+c,y:y+c]
q1 = q[x:x+c,y:y+c]
figure(1);imshow(z1, vmin=-vm, vmax=vm);title('Average')
figure(2);imshow(q1, vmin=-vm, vmax=vm);title('single')
print(z1.std());print(q1.std())


x =350;y = 300;c=100; vm=20e-9
figure(1);clf();imshow((qc[i])[x:x+c,y:y+c], vmin=-vm, vmax=vm)
figure(2);clf();imshow((qc[i+1])[x:x+c,y:y+c], vmin=-vm, vmax=vm)
figure(2);clf();imshow((qc[i+1])[x:x+c,y:y+c], vmin=-vm, vmax=vm)


# Moving/reading the OTT
#   - Reading the positions
parpos = ott.parabola.getPosition()
rmpos  = ott.referenceMirror.getPosition()
rmcar  = ott.referenceMirrorSlider.getPosition()
termp  = ott.temperature.getTemperature()

#   - Moving stuff
ott.referenceMirrorSlider.setPosition(975)

#  Calibration of the Alignment IntMat
#     - Command amplitudes and general initializations
par_piston = 0.7; 
par_tip = 100 ; par_tilt = 100
rm_tip  = 3   ; rm_tilt  = 3
command_amp_vector = np.array([par_piston, par_tip, par_tilt, rm_tip, rm_tilt])
nPushPull = 2; n_frames=20
#     - Acquistion
tnc = main.calibrate_PARAndRM(ott, interf, command_amp_vector, nPushPull, n_frames, delay=0)

tnc = '20230704_225757'
tnc = '20230714_133434'
tnc = '20230720_185546'
tnc = '20230720_192515'
# Alignment
#   - Definition of what and who
# alignment of focus
zern2corrf = np.array([2]);dofidf = np.array([0])# focus, parpist
#alignment of tilt onlu
zern2corrf = np.array([0,1]) ;dofidf = np.array([3,4])#TipTilt, RmTip, RmTilt
#alignment of tilt and coma
zern2corrf = np.array([0,1,3,4]); dofidf = np.array([1,2,3,4])# TipTilt, Coma,ParTip, ParTilt, RmTip, RmTilt
#zern2corrf = np.array([3,4]);dofidf = np.array([1,2,3,4])# Coma, ParTip, ParTilt, RmTip, RmTilt  !!wrong, not working

#   - Alignment
tna = main.align_PARAndRM(ott, interf, tnc, zern2corrf, dofidf, n_frames=1, doit=True)

# Measurements
#   - WFE and Zernike
img=interf.acquire_phasemap(5)
cc, m=zern.zernikeFit(img, [1,2,3,4,5,6,7,8,9,10])
print(cc[6:8])

#   - Time Series
tn = meas.opticalMonitoring(NMEAS,DELAY)
fl = th.fileList(tn)
imgave = th.averageFrames(0,100, fl)
imgave = th.removeZernike([1,2,3,4])

#   - Burst - Noise
phcamfocus.loadConfiguration(confnoise)
tau_vector = np.arange(1,100,2)
tn = interf.capture(2000)
interf.produce(tn)
dfpath = th.foldname.OPD_IMAGES_ROOT_FOLDER+'/'+tn+'/'
noise.convection_noise(dfpath, tau_vector)

#verification of Zernike computation during alignment
par0 = ott.parabola.getPosition()
par0[2]=par0[2]+0.2
ott.parabola.setPosition(par0)
time.sleep(5)
img=interf.acquire_phasemap(5)
cc, m=zern.zernikeFit(img, [1,2,3,4,5,6,7,8,9,10])
print(cc[3])
img=interf.acquire_phasemap(5)
cc, m=zern.zernikeFit(img, [1,2,3,4,5,6,7,8,9,10])
print(cc[3])
zern2corrf = np.array([2]) #focus
dofidf = np.array([0])# parpist
tna = main.align_PARAndRM(ott, interf, tnc, zern2corrf, dofidf, n_frames=4, doit=True)
img=interf.acquire_phasemap(5)
cc, m=zern.zernikeFit(img, [1,2,3,4,5,6,7,8,9,10])
print(cc[3])

#RMslider calibration
#limits: 0, upper position in 4D. 1328 lower position in 4D
rpos = ott.referenceMirrorSlider.getPosition()
ott.referenceMirrorSlider.setPosition(1000)

#  MEASUREMENTS
#Time Series
tn = meas.opticalMonitoring(NMEAS,DELAY)

#  Burst - Noise
phcamfocus.loadConfiguration(confnoise)
tau_vector = np.arange(1,100,2)
tn = interf.capture(2000)
interf.produce(tn)
dfpath = th.foldname.OPD_IMAGES_ROOT_FOLDER+'/'+tn+'/'
noise.convection_noise(dfpath, tau_vector)


#  Averaged frame for PAR subtraction tbc
from m4.misc import par_filterCenter as pf
tnott = '20230706_192127'  #same tracknum as in 20230707_test_par_registration
tnott = '20230715_113735' #with fans, plateu =17 nm
#aveimg, plist = pf.filterFrames(tnott,rad1=0.1,rad2=0.02,thr=1.4) #thr = 18 nm, 493 frames
fl = th.fileList(tnott)
aveimg = th.averageFrames(0,999, fl)
fname = th.foldname.OPD_SERIES_ROOT_FOLDER+'/'+tnott+ '/average.fits'
pyfits.writeto(fname, aveimg.data)
pyfits.append(fname, aveimg.mask.astype(np.uint8))



