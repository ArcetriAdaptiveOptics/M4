from m4.ground import read_data
from m4.mini_OTT import timehistory as th
from m4.ground import geo
from m4.utils import image_registration_lib as imgreg
from m4.utils.parabola_identification import ParabolaActivities
pa = ParabolaActivities()
#pa richiede un masked array ma lavora solo sulla maschera

tnm = '20240521_094351'
fname = '/mnt/backup/Archeopterix_20230517/Data/M4Data/OPTData/20240521_144221.fits'
q = read_data.readFits_data(fname)
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

zz1 = -1*zz+1
zm = -1*zz1*(-1*qm.mask+1)+1

qima = np.ma.masked_array(q,zm)
mpos = pa.rawMarkersPos(qima)
diam =5
npix = diam**2
athr = 3
mpos1 = pa.filterMarkersPos(mpos, npix/athr, npix*athr)

mpos1 = pa.filterMarkersPos(mpos, 10,100)


