import numpy as np
from photutils.centroids import centroid_2dg
from m4.ground import geo
from m4.utils import image_registration_lib as imgreg
center_act = 313

def findFrameCoord(imglist, actlist, actcoord):
    '''
    returns the position of given actuators from a list of frames
    '''
    pos = []
    for i in imglist:
        pos.append(findActuator(i))
    pos = (np.array(pos)).T

    frameCenter= imgreg.marker_general_remap(actcoord[:,actlist],
                                             pos,
                                             actcoord[:,(center_act,center_act)])
    #the last variable has been vectorized (by adding a second element) don't know why but so it works
    frameCenter = frameCenter[:,0]
    return frameCenter

def findActuator(img):
    '''
    Finds the coordinates of an actuator, given the image with the InfFunction masked around the act.
    img: masked array
        image where the act is to be searched
    Return
    imgout: array
        coordinates of the act
    '''
    imgw    = extractPeak(img,radius=50)
    pos     = centroid_2dg(imgw)
    return pos



def extractPeak(img, radius = 50):
    '''
    Extract a circular area around the peak in the image
    '''
    yp, xp  = np.where(img == np.max(abs(img)))
    img1    = img*np.invert(img.mask)
    m       = np.invert(geo.draw_mask(img.mask, yp, xp, radius))
    imgout  = np.ma.masked_array(img1, m)
    return imgout

def combineMasks(imglist): #!!! deve sparire
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
    mmask = product(mm, 0)
    return mmask


#how to find the M4 segment position in frame:
# Sc: segment center, coordinate of the central actuator in segment (defined)
# Scf: nominal position in the frame of the Sc posistion (this is defined by moving the truss to aim a given segment
# procedure: we acquire a number of actuators and measure their position in the frame (MAct);
# then we build the system with nominal coordinate (PAct) and measured ones (MAct).
# the current Scf_i is found with:
    #   Scf_i = imgreg.marker_general_remap(MAct, PAct, Scf) (check the inputs)
