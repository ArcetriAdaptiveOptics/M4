from m4.ground import read_data
from m4.mini_OTT import timehistory as th
from m4.ground import geo
from m4.utils import image_registration_lib as imgreg
from m4.utils.parabola_identification import ParabolaActivities
pa = ParabolaActivities()
#pa richiede un masked array ma lavora solo sulla maschera

tnm = '20240521_094351'
fl=th.fileList(tnm)
aveimg = th.averageFrames(0,2,fl)
fname = '/mnt/backup/Archeopterix_20230517/Data/M4Data/OPTData/20240521_144221.fits'
q = read_data.readFits_data(fname)
#here I need to create a test image
qmask = np.ones(q.shape)
qmask[q <2*q.mean()]=0
q = np.ma.masked_array(q, -1*qmask+1)

#now i need to create the test mask-inverted markers
pos = imgreg.getMarkers(tnm, diam=28)
rad=30
mm = np.ones([2000,2000])
zz = mm
for i in range(0,9):
    zz = zz*geo.draw_mask(mm,pos[0,i],pos[1,i],rad, out=1)

m1 = zz *q.mask
zz1 = -1*zz+1
q1m = -1*q.mask+1
m1 = zz1*q1m
imshow(m1)

qimg = np.ma.masked_array(q.data,-1*m1+1)  #this is our test image

#testing the markerd identification
mpos = pa.rawMarkersPos(qimg)
diam =5
npix = diam**2
athr = 3
mpos1 = pa.filterMarkersPos(mpos, npix/athr, npix*athr)

#finding the circles
mlist0 = [2,3,4]
mlist1 = [0,1,5,8]
c0, axs0, r0 = pa._fitEllipse(mpos1[1,mlist0],mpos1[0,mlist0]); c0 = np.real(c0)
c1, axs1, r1 = pa._fitEllipse(mpos1[1,mlist1],mpos1[0,mlist1]); c1 = np.real(c1)





markmask = -1*aveimg.mask+1
imgmask = aveimg.mask
dmask = qmask*imgmask



q1 = q.copy()
thr = 30
rad = 50
q[q < thr]=0
imshow(q)
qm = np.ma.masked_array(q,q==0)
pos = imgreg.getMarkers(tnm, diam=28)

mm = np.ones([2000,2000])
zz = mm
for i in range(0,9):
    zz = zz*geo.draw_mask(mm,pos[0,i],pos[1,i],rad, out=1)

#qima = np.ma.masked_array(q,zz) #this is wrong. we need the detector mask to be masked in the hole!!

zz1 = -1*zz+1
zm = -1*zz1*(-1*qm.mask+1)+1

qima = np.ma.masked_array(q,zm)
mpos = pa.rawMarkersPos(qima)
diam =5
npix = diam**2
athr = 3
mpos1 = pa.filterMarkersPos(mpos, npix/athr, npix*athr)

mpos1 = pa.filterMarkersPos(mpos, 10,100)


