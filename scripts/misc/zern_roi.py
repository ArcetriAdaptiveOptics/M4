
import numpy as np
from opticalib.ground import geo
from opticalib.ground import zernike as zern
from opticalib.ground.osutils import load_fits
from opticalib.ground import roi

## this part is for initializing a test session ###
fname = '/home/labot/mode_0114.fits'
img = load_fits(fname)
imgroi = roi.roiGenerator(img, 2)
mm0 = np.invert(img.mask)
mm = geo.draw_mask(img.data*0,460,460,450,out=0)
imshow(mm+img)
#here we suppose that zernikeAuxmask works


def roizern(img, z2fit, auxmask =None, roiid=None, local =True):
    '''
    to compute Zernike over all the passed roiid.
    local  --> ALL the zernike (computed over all the rois) are returned
    global --> the average of the Zernike is returned
    '''
    if roiid is not None:  #
        roiimg = roi.roiGenerator(img) #non è disponibile un parametro passato per dire QUANTE roi cercare. funziona anche senza?
        nroi = len(roiid)
    else:
        nroi=1
    zcoeff = np.zeros([nroi, len(z2fit)])
    zsurf  = []
    for i in range(nroi):
        cc, _ =zern.zernikeFitAuxmask(roiimg[i], auxmask=auxmask, z2fit)
        zcoeff[i,:] = cc
    if local is False:
        zcoeff = zcoeff.mean(axis=0)

    return zcoeff

def tiltDetrend(img, auxmask, roi2Calc, roi2Remove):
    '''
    computes the Zernikes (PTT only)  over the roi2Calc, then produces the corresponding shape over the roi2Remove mask, and subtract it.
    '''
    roiimg = roi.roiGenerator(img)
    zcoeff = roizern(img, [1,2,3], auxmask, roiid=roi2Calc, local=False) #returns the global PTT evaluated over the roi2Calc areas
    _, zmat = zern.zernikeFit(auxmask,[1,2,3]) #returns the ZernMat created over the entire circular pupil
    surf2Remove = zern.zernikeSurface( auxmask,zcoeff, zmat)
    surf2Remove = np.ma.masked_array(surf2remove.data, roi2Remove.mask)
    detrendedImg = roi2remove - surf2Remove
    return detrendedImg


    




