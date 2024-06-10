from photutils.centroids import centroid_2dg
from m4.ground import geo


from astropy.io import fits as pyfits

a= '/home/runa/act-marker.fits'

hd=pyfits.open(a)
q = hd[0].data
img = np.ma.masked_array(q[0,:,:],-1*(q[1,:,:])+1)

nmask = -1*img.mask+1
yp, xp = np.where(img == np.max(abs(img)))
img1 = geo.draw_mask(nmask,yp,xp, 50,out=1)

newmask = nmask * img1
xx = centroid_2dg(img)

