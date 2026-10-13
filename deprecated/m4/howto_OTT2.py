#in git/labott: calpyott 

import m4
from m4.devices import opt_beam
from opticalib.core import root as fname
ott,_,interf=m4.create_ott()

par = opt_beam.Parabola(ott)
refmirr = opt_beam.ReferenceMirror(ott)

from m4 import userscripts
myott = userscripts.OTTScripts(ott, interf)

myott.config4D4Alignment()
myott.alignOTT(1,0)
myott.acquireTimeAverage(3,1)

myott.configureOTTrefMirrorOut()
myott.configureOTT4Segment()
myott.whereisRef()
myott.whereisBeam()


from scripts.misc import ott_measurements as mm
meas=mm.Measurements(interf, ott)
tn=meas.opticalMonitoring(3,0)

al =  m4.OttAligner(ott, interf)
cmdAmp = [0.7, 100, 100, 6,6]; n_frames =5 #PAR + RM
cmdAmp = [0,0,0,0,0,0.01,0.01]; n_frames=5  #DP only, 7x7
al.calibrate_alignment(cmdAmp, n_frames,save=True)
tncal ='20250923_142114'  #tncal='20250909_141737'  #PAR+RM only
tncal = '20250924_151216'  DP only, DoF [5,6]
al.load_calibration(tncal)
al.correct_alignment([3,4],[0,1],applycmd=False, n_frames=1) #DoF, Zern

rm= opt_beam.ReferenceMirror(ott)
rm?
rm.rmTipTilt([1,1])   #diff tilt command to RM

#analyze frames
#fl = osu.getFileList(tns[0], fold='INTMatrices', key='Cube')
#in OPDSeries, maps
f = osu.getFileList(tn, key='20')

#from micws 
#module load dp
## da git/M4
##    calpy -f /data/Arcetri/Data/M4Data/

import numpy as np
import m4
from m4.devices import opt_beam
from opticalib.core import root as fname
ott,_,interf=m4.create_ott()

from m4 import userscripts
myott = userscripts.OTTScripts(ott, interf)

from opticalib import PhaseCam, AdOpticaDm
from opticalib.dmutils import iff_module as ifm, iff_processing as ifp, iff_acquisition_preparation as ifa
from opticalib.dmutils.flattening import Flattening
from opticalib import load_fits, save_fits
from matplotlib.pyplot import *
ion()
%matplotlib qt

dm = AdOpticaDm()

interf = PhaseCam(6110)

c3=load_fits('/data/Arcetri/Data/M4Data/SysConfig/20250717_150000_flatId0.fits')
c3 = load_fits('/data/Arcetri/Data/M4Data/SysConfig/20250909_155400.fits')  #flat 30 modes and local tilt
c4 = load_fits('/data/Arcetri/Data/M4Data/SysConfig/20250930_155500_pos.fits')
c5 = load_fits('/data/Arcetri/Data/M4Data/SysConfig/20250930_172600_flat.fits')
dm.set_shape(c3*0.1) #....

#ifc = ifa.IFFCapturePreparation(dm)
#thist =ifc.createTimedCmdHistory(np.arange(10)+111,200e-9)


tn=ifm.iffDataAcquisition(dm, interf, [0,1,2,3,4,5], 1e-6)  #cmdOffset
interf.produce(tn)

ifp.process(tn)# !!!! when nworking multi-segment, the individual piston may produce a false Trigger Found event. we need to implement a local pist-tiptilt fitting
cube = ifp.saveCube(tn)

f = Flattening(tn)
f.filterIntCube([1,2,3])

img = interf.acquire_map()
f.loadImage2Shape(img)
f.computeRecMat(3)  #nmodes 2 remove
fcmd = f.computeFlatCmd(5)

# test of eigenmodes inspection
from opticalib.ground import computerec as _crec
rec = _crec.ComputeReconstructor

###log 20250910

surface = al._surface.copy()
from opticalib.ground import geo
from opticalib.ground import zernike as zern
from opticalib.ground import roi
surface = al._surface.copy()
cir = geo.qpupil(-1*surface.mask + 1)
mm=geo.draw_mask(surface.data * 0, cir[0], cir[1], 1.44 / 0.00076 / 2, out=0)
img = al._acquire[0](nframes = 1)
tt=roizern(img,[1,2,3],mm,[0,1])

#per processare eliminando una roi:
ifp.process(tn, roi=0)    #roi=0 corresponds to trigger mode on shell at ROI0


#flattening automated
tn='20250911_062742'
f=dmutils.flattening(tn)
f.applyFlatCommand(dm, interf, modes2flat=222,modes2discard=1, cmdOffset=cmd, nframes=5)


def roizern(img, z2fit, auxmask =None, roiid=None, local =True):
    if roiid is not None:  #
        roiimg = roi.roiGenerator(img) #non è disponibile un parametro passato per dire QUANTE roi cercare. funziona anche senza?
        nroi = len(roiid)
    else:
        nroi=1
    if auxmask is None:
        auxmask2use = img.mask
    else:
        auxmask2use = auxmask
    zcoeff = np.zeros([nroi, len(z2fit)])
    zsurf  = []
    for i in range(nroi):
        img2fit = np.ma.masked_array(img.data, roiimg[i])
        cc, _ =zern.zernikeFitAuxmask(img2fit, auxmask2use, z2fit)
        zcoeff[i,:] = cc
    if local is False:
        zcoeff = zcoeff.mean(axis=0)
    return zcoeff

c3 = load_fits('/data/Arcetri/Data/M4Data/SysConfig/20250909_155400.fits')  #flat 30 modes and local tilt
dm.set_shape(c3*0.1)

 
