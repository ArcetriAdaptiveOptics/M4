import numpy as np
from m4.configuration import start
from m4.ground import zernike as zern
from m4.misc import par_meas as pp
from m4.ground import timestamp
import time
from m4 import main
from m4.mini_OTT.measurements import Measurements

conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
ott, interf, dm = start.create_ott(conf)

ts = timestamp.Timestamp()

tnc = '20230113_102942'

zern2corrf = np.array([0,1,2]) #TipTilt focus
dofidf = np.array([0,1,2])# parpist, ParTip, ParTilt
par0=ott.parabola.getPosition()

tna = main.align_PARAndRM(ott, interf, tnc, zern2corrf, dofidf,n_frames=4)

meas = Measurements(ott, interf)
delay=60
nmeas = 1000
tn = meas.opticalMonitoring(nmeas, 60)



#to acquire a PAR image and compute zern
img1,tn = pp.acq_par(16,1)
cc = pp.get_zern(img1)
#to evaluate the coma
img,tn= pp.acq_par(16,1);print(tn);cc = pp.get_zern(img);print(cc[6:8])

img,tn= pp.acq_par(16,1);print(tn);cc = pp.get_zern(img);print('Astigm');print(cc[4:6]);print('Coma');print(cc[6:8])



#calibrating the alignment intmat
par_piston = '?'
par_tip = '?'
par_tilt = '?'
rm_tip = '?'
rm_tilt = '?'
command_amp_vector = np.array([par_piston, par_tip, par_tilt, rm_tip, rm_tilt])
main.showCommandMatrixBeforeCalibration(command_amp_vector)

nPushPull = 1; n_frames=9
tnc = main.calibrate_PARAndRM(ott, interf, command_amp_vector, nPushPull, n_frames, delay=0)

# id for tip tilt
zern2corr = np.array([0,1]) #TipTilt
dofid = np.array([1,2])# ParTip, ParTilt

zern2corrf = np.array([0,1,2]) #TipTilt focus
dofidf = np.array([0,1,2])# parpist, ParTip, ParTilt

main.showCommandForParAndRmBeforeAlignement(ott, interf, tnc, 2, zern2corr, dofid)

tnc = '20230113_102942'
tna = main.align_PARAndRM(ott, interf, tnc, zern2corrf dofidf,n_frames=4)

#test of shaking the PAR
dp = np.array([0,0,1,0,0,0])
par0 = ott.parabola.getPosition()
ott.parabola.setPosition(par0+dp)
ott.parabola.setPosition(par0-dp)
tnc = main.calibrate_PARAndRM(ott, interf, command_amp_vector, nPushPull, n_frames, delay=0)
img,tn= pp.acq_par(16,1);print(tn);cc = pp.get_zern(img);print('Astigm');print(cc[4:6]);print('Coma');print(cc[6:8])


#burst
dd = 'D:/M4/Capture/'
tn = ts.now()
interf.burstFramesToSpecificDirectory(dd+tn+'/', 1000)


