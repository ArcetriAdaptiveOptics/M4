%run /home/m4/git/M4/m4/misc/initOTT.py
#  4D Configurations
confcenter = 'D:/config/20230714_ott_center.4Dini'
confnoisecenter = 'D:/config/20230714_ott_center-LowRes.4Dini'
confnomask = 'D:/config/20230706_ott_automask.4Dini'
confrm1328 = 'D:/config/20230707_ott_RM1328.4Dini'
confnomasklowres = 'D:/config/20230706_ott_automask-lowRes.4Dini'
confnoise = 'D:/config/20230706_ott_center-lowResNoise.4Dini'
confnoisezoom = 'D:/config/20230706_ott_center-lowRes-Zoom-Noise.4Dini'

phcamfocus.loadConfiguration(confnomask)

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
rm_tip  = 6   ; rm_tilt  = 6
command_amp_vector = np.array([par_piston, par_tip, par_tilt, rm_tip, rm_tilt])
nPushPull = 2; n_frames=20
#     - Acquistion
tnc = main.calibrate_PARAndRM(ott, interf, command_amp_vector, nPushPull, n_frames, delay=0)

tnc = '20230704_225757'
tnc = '20230714_133434'
tnc = '20230720_185546'
tnc = '20230720_192515'
tnc = '20231027_095920' #!! with global modes, 20 frames, large amplitudes
tnc = '20231027_131458' #global modes, high SNR
# Alignment
#   - Definition of what and who
# alignment of focus
zern2corrf = np.array([2])       ;dofidf = np.array([0])# focus, parpist
#alignment of tilt onlu
zern2corrf = np.array([0,1])     ;dofidf = np.array([3,4])#TipTilt, RmTip, RmTilt
#alignment of tilt and coma
zern2corrf = np.array([0,1,3,4]) ; dofidf = np.array([1,2,3,4])# TipTilt, Coma,ParTip, ParTilt, RmTip, RmTilt
#zern2corrf = np.array([3,4]);dofidf = np.array([1,2,3,4])# Coma, ParTip, ParTilt, RmTip, RmTilt  !!wrong, not working

#   - Alignment
tna = main.align_PARAndRM(ott, interf, tnc, zern2corrf, dofidf, n_frames=1, doit=True)

# Measurements
#   - WFE and Zernike
img=interf.acquire_phasemap(5)
cc, m=zern.zernikeFit(img, [1,2,3,4,5,6,7,8,9,10])
print(cc[6:8])

#   - Time Series
NMEAS = 1
DELAY = 1
tn = meas.opticalMonitoring(NMEAS,DELAY)
fl = th.fileList(tn)
imgave = th.averageFrames(0,100, fl)
imgave = th.removeZernike([1,2,3,4])

#   - Burst - Noise
phcamfocus.loadConfiguration(confnoise)
tn = interf.capture(2000)
interf.produce(tn)
dfpath = th.foldname.OPD_IMAGES_ROOT_FOLDER+'/'+tn+'/'
tau_vector = np.arange(1,100,2)
noise.convection_noise(dfpath, tau_vector)
template=np.array([3,11,25,37,51])
noise.noise_vibrations(dfpath,template)

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
# burst analysis, results saved in "/mnt/m4storage/Data/M4Data/Results/Noise"
dfpath = os.path.join(th.foldname.OPD_IMAGES_ROOT_FOLDER,tn)
noise.convection_noise(dfpath, tau_vector)
template=np.array([3,11,25,37,51])
noise.noise_vibrations(dfpath,template)


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




