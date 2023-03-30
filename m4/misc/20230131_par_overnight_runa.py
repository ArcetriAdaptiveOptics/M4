from astropy.io import fits as pyfits
from m4.mini_OTT import timehistory as th
from matplotlib import pyplot as plt
import numpy as np

tn='20230127_163227'
tn='20230129_100739'  # day!!
tn='20230131_181616'  # cover
tn='20230201_183722'  # cover
tn='20230202_175230'  # cover, fan
tn='20230203_182356'  # cover, fan, box on
tn='20230204_160455'  # cover, fan, box on


fl = th.fileList(tn)
nf = len(fl)
tvec = th.timevec(tn)
tvec=(tvec-tvec[0])*24
base = '/mnt/m4storage/Data/M4Data/OPTData/OPD_series/'
zname=base+tn+'/zernike.fits'
z=pyfits.open(zname)[0].data
a5=z[:,5]
a6 = z[:,6]
plt.plot(th.runningMean(a5,10),'.'); plt.title(tn+' Ast-Runnmean')



q0=th.averageFrames(100,800,fl)
i=0
plt.imshow(th.removeZernike(th.frame(i, fl)-q0,[1,2,3,4]))
rr = []
for i in range(nf):
    qi = th.removeZernike(th.frame(i,fl)-q0,[1,2,3,4,7,8])
    rr.append(qi.std())

plt.plot(th.runningMean(rr,10), '.'); plt.title(tn+' RMS-ref')

zz = th.zernikePlot(fl)
nn = 100
dd = th.runningDiff(tn,2)
plt.plot(tvec[4::2], dd, '.');plt.title(tn+': RunningDiff-RMS,2min'); plt.xlabel('Time since start [hr]')
q=[]
for i in range(10):
    qq=th.removeZernike(th.averageFrames(i*100,i*100+nn-1,fl))
    q.append(qq)

plt.clf();plt.imshow(q[1]-q[0]);plt.colorbar();plt.title(tn+', avg diff I'); np.std(q[1]-q[0])
plt.clf();plt.imshow(q[4]-q[3]);plt.colorbar();plt.title(tn+', avg diff II'); np.std(q[4]-q[3])
plt.clf();plt.imshow(q[8]-q[7]);plt.colorbar();plt.title(tn+', avg diff III'); np.std(q[8]-q[7])

nn = 100
step=int(nf/nn)
q=[]
for i in range (nn):
    q.append(th.removeZernike(th.frame(i*step,fl)))
qa=np.ma.masked_array(q)
stimg = np.std(qa,0)
plt.clf();plt.imshow(stimg, vmax=50e-9);plt.colorbar();plt.title(tn+' PixStddev')

# lab =['ciao', ciao']
#how to draw labels: plot(ff,vv[:,i], label=lab[i])



