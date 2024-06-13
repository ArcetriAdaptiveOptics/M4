from photutils.centroids import centroid_2dg
from m4.ground import geo


from astropy.io import fits as pyfits

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
