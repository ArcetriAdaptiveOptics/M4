from scripts.misc.IFFPackage import actuator_identification_lib as fa
from m4.utils import image_registration_lib as imgreg
from photutils.centroids import centroid_2dg
from m4.ground import geo
from m4.ground import read_data as rr
from m4.devices import deformable_mirror as dm
from astropy.io import fits as pyfits
m4u=dm.M4AU()

fname = '/mnt/cargo/data/M4/Data/M4Data/SimOPTData/OPDImages/20160516_114916/img_0008.fits'
q = rr.read_phasemap(fname)

pos = fa.findActuator(img)

b = '/mnt/cargo/data/M4/Data/M4Data/SimOPTData/OPDImages/20160516_114916/img_000'
f = ['7','8','9']
imglist=[]
for i in f:
    imglist.append(rr.read_phasemap((b+i+'.fits'))


actcoord = m4u.actCoord
actlist = [ 585,110,609]
fa.findFrameCoord(imglist, actlist, actcoord)


#------

a= '/home/runa/act-marker.fits'

hd=pyfits.open(a)
q = hd[0].data
img = np.ma.masked_array(q[0,:,:],-1*(q[1,:,:])+1)

nmask = -1*img.mask+1
yp, xp = np.where(img == np.max(abs(img)))
m1 = geo.draw_mask(nmask*0,yp,xp, 50)
newmask = nmask*m1
img1 = np.ma.masked_array(img,-1*newmask+1)

xx = centroid_2dg(img1)




def combineMasks(imglist):
    '''
    combine masks layers of masked arrays, or a list of masks, to produce the intersection masks: not masked here AND not mnaked there
    masks are expected as in the np.ma convention: True when not masked
    return:
        intersection mask

    '''
    imglistIsMaskedArray = True
    imglistIsmasksList   = False
    mm = []
    for i in imglist:
        if imglistIsMaskedArray:
            mm.append(np.invert(i.mask).astype(int))
        if imglistIsmasksList:
            mm.append(np.invert(i).astype(int))


    mmask = product(mm,0)
    return mmask
