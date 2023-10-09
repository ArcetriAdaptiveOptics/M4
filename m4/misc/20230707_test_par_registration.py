
import numpy as np
from matplotlib import pyplot as plt
from m4.mini_OTT import timehistory as th
from m4.utils.parabola_identification import ParabolaActivities
from matplotlib.patches import Circle

'''
OTT
TN
20230707_142156
#markers visible in OTT center view:
    np.array([6,7,8,11,12,13,14,17,18])

CGH

2000x2000, offs 28,0
TN          	DIM MARKERS [MM]	N_MEAS
20230328_175944	20	                20
20230328_180252	20	                20
20230328_181123	10	                20
20230328_181333	10	                20

1680x1680, offs 180,180
TN	            DIM MARKERS [MM]	N_MEAS
20230331_181519	20	                20
20230331_181658	20	                20
20230331_181946	10	                20
20230331_182118	10	                20

relative offset is 180-28, 180
'''
pa = ParabolaActivities()
offs = [180-28,180]
nmark = 25
offvec = np.zeros([2,nmark])
for i in range(nmark):
    offvec[:,i]=offs

circ0 = [8,11,12]
circ1 = [6,7,13,14,17,18]
circ2 = [3,4,5,15,16,20,21]
circ3 = [0,1,2,9,10,19,22,23,24]
radii = np.array([0.1, 0.2505, 0.452, 0.6911])

def coord2ima(cc):
    cc1=np.array([cc[1,:],cc[0,:]])
    #cc2=np.array([2000-cc[1,:],cc[0,:]])
    #cc3=np.array([cc[1,:],2000-cc[0,:]])
    #cc4=np.array([2000-cc[1,:],2000-cc[0,:]])

    return cc1

def getMarkers(tn,flip=False, diam=24, thr=0.2):
    npix = 3.14*(diam/2)**2
    fl = th.fileList(tn)
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
# OTT registration, RefMirror vs CGH
tn0 = '20230328_175944' #2000x2000
tn1 = '20230707_142156'
fl0 = th.fileList(tn0)
img0=th.frame(0, fl0)
img0= np.fliplr(img0)
fl1 = th.fileList(tn1)
img1=th.frame(0, fl1)

p0 = getMarkers(tn0,flip=True, diam=24)
p1 = getMarkers(tn1, diam=28)
p0 = coord2ima(p0)
p1 = coord2ima(p1)
plt.imshow(img0)
for i in range(np.shape(p0)[1]):
    plt.text(p0[0,i],p0[1,i],i)
plt.figure(2)
plt.imshow(img1)
for i in range(np.shape(p1)[1]):
    plt.text(p1[0,i],p1[1,i],i)

#CGH vs OTT markers match, at Center_view
mark_cgh = np.array([8,11,12,6,13,18,17,14,7])
mark_ottc= np.array([2,3,4,1,6,8,7,5,0])

pcgh = p0[:,mark_cgh]
pott = p1[:,mark_ottc]


 
