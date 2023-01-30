'''
Authors
  - ?. ?: written in 2022
'''

#da far girare prima di incominciare le misure
import numpy as np
from matplotlib import pyplot as plt
from m4.misc import par_meas as pp
from m4.ground import read_data as rr
from m4.mini_OTT import timehistory as th 
from m4.mini_OTT.measurements import Measurements 
from m4 import main
from astropy.io import fits as pyfits
from m4.configuration import start
conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
ott, interf, dm = start.create_ott(conf)
meas = Measurements(ott,interf)

#acquisition of time series
NMEAS = 200
DELAY = 5 
tn = meas.opticalMonitoring(NMEAS,DELAY)

tnc = '20230113_102942'
zern2corrf = np.array([0,1,2]) #TipTilt focus
dofidf = np.array([0,1,2])# parpist, ParTip, ParTilt

tna = main.align_PARAndRM(ott, interf, tnc, zern2corr, dofid, n_frames=4)


#acquisition of single (averaged) frame
img1,tn = pp.acq_par(16,1)
cc = pp.get_zern(tn)

fname='/mnt/m4storage/Data/M4Data/OPTData/PARTest/20221220/20221220_125418.fits'
fname='/mnt/m4storage/Data/M4Data/OPTData/PARTest/20221220/20221220_101045.fits'
sname = '20221220_143412'

fold = th.findTracknum(sname)
flist = th.fileList(sname)
zz = th.zernikePlot(flist[0:10])

img = th.read_phasemap(fname)

#analysis of time series
w   = '/mnt/m4storage/Data/M4Data/OPTData/OPD_series/'
tn = '20230113_133110'
tn = '20230113_143700'
zz = (pyfits.open(w+tn+'/zernike.fits'))[0].data
zm = np.mean(zz,0)
zs = np.std(zz,0)
print(zm[5:7]);print(zs[5:7])
 #20230113_133110
#average:     [-2.5e-08 -4.6e-09]
#dispersion:  [4.1e-09 4.0e-09]
#20230113_143700
#average:    [-2.4e-08 -3.6e-09]
#dispersion: [3.8e-09 3.5e-09]
n2 = int(len(zz)/2)
m0 = (np.mean(zz[0:n2,5]),np.mean(zz[n2:,5]))
m0 = (np.mean(zz[0:n2,6]),np.mean(zz[n2:,6]))


#analysis of time series
def meas_ast(tn):
    fl=th.fileList(tn)
    q=th.averageFrames(0,99,fl)
    c,m=th.zernike.zernikeFit(q,[1,2,3,4,5,6,7,8,9,10])
    return c[4:6]





#PAR meas Zernike.fits file analysis
from astropy.io import fits as pyfits
#import glob
where = '/mnt/m4storage/Data/M4Data/OPTData/OPD_series/'
fold = '20221220_143412'
fname = 'zernike.fits'
file_name = where + fold + '/' + fname
hdul = pyfits.open(file_name)
data = hdul[0].data
time = data[:,0]
coeff = data[:,1::]
ave = np.average(coeff,0)
sstd = np.std(coeff,0)
plt.figure(figsize=(8,5))
plt.plot(coeff[:,4],'-o', label='ast X')
plt.plot(coeff[:,5],'-o', label='ast Y')
plt.legend(); plt.grid(); plt.ylabel('nm RMS'); plt.xlabel('frame #'); plt.title(fold + ' Astigmatism plot')
plt.figure(figsize=(8,5))
plt.plot(coeff[:,6],'-o', label='coma X')
plt.plot(coeff[:,7],'-o', label='coma Y')
plt.legend(); plt.grid(); plt.ylabel('nm RMS'); plt.xlabel('frame #'); plt.title(fold + ' Coma plot')
print('Average Astigmatism =  %.3e ±  %.3e, %.3e ± %.3e  m' %(ave[4], sstd[4], ave[5], sstd[5]))
print('Average Coma =  %.3e ± %.3e, %.3e ± %.3e m' %(ave[6],sstd[6], ave[7], sstd[7]))




import glob
path = '/mnt/m4storage/Data/M4Data/OPTData/PARTest/'
fold = '20221220'
where = path + fold + '/'
D = sorted(glob.glob(where + '*.fits'))
coeff = []
for ii in range(len(D)):
    img = th.read_phasemap(D[ii])
    cc = pp.get_zern(img)
    coeff.append(cc)
coeff = np.array(coeff)





#D = sorted(glob.glob(where + '20221220_*'))
#coeff_series = []
#for ii in range(len(D)):
#    img = th.read_phasemap(D[ii])
#    cc = pp.get_zern(img)
#    coeff_series.append(cc)

#analysis of time histories
tnlist = ['20230125_131542','20230125_144054','20230125_152344','20230125_171814']
zz = []
img = []
dd =[]
dq=[]
for i in tnlist:
    fl=th.fileList(i)
    nf=len(fl)
    q= th.averageFrames(0, nf, fl)
    img.append(q)
    zf, m= th.zernike.zernikeFit(q, [1,2,3,4,5,6])
    zz.append(zf)
    d1=th.runningDiff(i,1)
    dd.append(d1)
    dq0=th.averageFrames(0,int(nf/2),fl)-th.averageFrames(int(nf/2),nf,fl)
    dq0=th.removeZernike(dq0)
    dq.append(dq0)
    #for jj in np.arange(nf-50):
    #    th.averagedFrames(jj,jj+50,fl)


tn = pp.capture(1000)
pp.produce(tn)


