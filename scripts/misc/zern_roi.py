
import numpy as np
from opticalib.ground import geo
from opticalib.ground import zernike as zern
from opticalib.ground.osutils import load_fits
from opticalib.ground import roi

from arte.types.mask import CircularMask
from arte.utils.zernike_generator import ZernikeGenerator


## this part is for initializing a test session ###
fname = '/home/labot/mode_0114.fits'
img = load_fits(fname)
ss = shape(img)
imgroi = roi.roiGenerator(img, 2)
mm0 = np.invert(img.mask)
#mm = geo.draw_mask(img.data*0,460,460,450,out=0)

img = np.ma.masked_array(img.data, imgroi[0])

mm=geo.draw_mask(np.zeros([920,920]),460,460,450,out=0)
imshow(mm+img)
cc, mat =zern.zernikeFitAuxmask(img, mm, [1,2,3])

#creating a differential tilt
mmask = np.ma.masked_array(mm, mask=mm==0)
cmask = CircularMask.fromMaskedArray(mmask, mask=mm==0)
zgen = ZernikeGenerator(cmask)
tip = zgen.getZernike(2)
tilt= zgen.getZernike(3)
coeff0 = [1,4,6]
coeff1 = [5,10,2]
img0 = np.ma.masked_array(tip*coeff0[1]+tilt*coeff0[2]+coeff0[0], imgroi[0])
img1 = np.ma.masked_array(tip*coeff1[1]+tilt*coeff1[2]+coeff1[0], imgroi[1])
allmask = imgroi[0].copy()
allmask[imgroi[1]==0]=0
imgf = np.zeros(ss)
imgf[imgroi[0]==0]=img0.data[imgroi[0]==0]
imgf[imgroi[1]==0]=img1.data[imgroi[1]==0]
imgf = np.ma.masked_array(imgf, allmask)
imshow(imgf)



mmask = np.ma.masked_array(mm, mask=mm==0)
cmask = CircularMask.fromMaskedArray(mmask, mask=mm==0)
zgen = ZernikeGenerator(cmask)
w = zgen.getZernike(2)
imshow(w)

cc, mat =zern.zernikeFitAuxmask(w, mm, [1,2,3])
#here we suppose that zernikeAuxmask works

q=roizern(imgf, [1,2,3],mm, [0,1])

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

v=tiltDetrend(imgf,mm, [0],[1])

def tiltDetrend(img, auxmask, roi2Calc, roi2Remove):
    '''
    computes the Zernikes (PTT only)  over the roi2Calc, then produces the corresponding shape over the roi2Remove mask, and subtract it.
    '''
    roiimg = roi.roiGenerator(img)
    zcoeff = roizern(img, [1,2,3], auxmask, roiid=roi2Calc, local=True) #returns the global PTT evaluated over the roi2Calc areas
    am = np.ma.masked_array(auxmask, auxmask==0)
    _, zmat = zern.zernikeFit(am,[1,2,3]) #returns the ZernMat created over the entire circular pupil
    surf2Remove = zern.zernikeSurface( am,zcoeff[roi2Calc,:], zmat)
    #surf2Remove[roiimg[roi2Remove==0]] =0
    #surf2Remove = np.ma.masked_array(surf2remove.data, roi2Remove.mask)
    detrendedImg = imgf - surf2Remove
    return detrendedImg


    




