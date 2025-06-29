# -*- coding: utf-8 -*-
"""
Autors
  - R. Briguglio: created on Mon Mar 16 11:00:08 2020
  - F Miceli: add functionality on march 2022

"""
import numpy as np
from matplotlib import pyplot as plt
from scipy import ndimage
from skimage.measure import EllipseModel
#from skimage.draw import ellipse
from skimage.measure import CircleModel
from skimage.draw import disk
import scipy.interpolate


def qpupil_circle(image, pixel_dir=0):
    '''
    Function for...
    Created by Federico
    NOTA: la funzione usa come standard la direzione y per determinare la dimensione dei pixel

    pixel_dir: int
        indicates which direction to use for counting the number of pixels in the image.
        Y direction as standard
    '''
    aa = np.shape(image)
    imagePixels = aa[pixel_dir] #standard dir y
    circ = CircleModel()
    cnt = _trova_punti_bordi_ima2(image, imagePixels)
    circ.estimate(cnt)
    xc, yc, radius = np.array(circ.params, dtype=int)
    maskedd = np.zeros((imagePixels, imagePixels), dtype=np.uint8)
    rr, cc = disk((xc, yc), int(radius))
    maskedd[rr, cc] = 1

    idx = np.where(maskedd==1)
    ss = np.shape(maskedd)
    x = np.arange(ss[0]).astype(float)
    x = np.transpose(np.tile(x, [ss[1], 1]))
    y = np.arange(ss[1]).astype(float)
    y = np.tile(y, [ss[0], 1])
    xx = x
    yy = y
    #maxv = max(xx[idx])
    #minv = min(xx[idx])
    xx = xx - xc
    xx = xx/radius
    #maxv = max(yy[idx])
    #minv = min(yy[idx])
    yy = yy - yc
    yy = yy/radius

    return xc, yc, radius, xx, yy

def qpupil_ellipse(image, pixel_dir=0):
    '''
    Function for...
    Created by Federico
    NOTA: la funzione usa come standard la direzione y per determinare la dimensione dei pixel
    '''
    aa = np.shape(image)
    imagePixels = aa[pixel_dir] #standard dir y
    ell = EllipseModel()
    cnt = _trova_punti_bordi_ima2(image, imagePixels)
    ell.estimate(cnt)
    xc, yc, a, b, theta = np.array(ell.params, dtyp =int)
    maskedd = np.zeros((imagePixels, imagePixels), dtype=np.uint8)
    radius = max(a, b)
    rr, cc = disk((xc, yc), int(radius))
    maskedd[rr, cc] = 1

    idx = np.where(maskedd==1)
    ss = np.shape(maskedd)
    x = np.arange(ss[0]).astype(float)
    x = np.transpose(np.tile(x, [ss[1], 1]))
    y = np.arange(ss[1]).astype(float)
    y = np.tile(y, [ss[0], 1])
    xx = x
    yy = y

    xx = xx - xc
    xx = xx/radius

    yy = yy - yc
    yy = yy/radius

    return xc, yc, radius, xx, yy

def _trova_punti_bordi_ima2(image, imagePixels):
    '''
    Function for...
    Created by Federico
    '''
    x = image
    val = []

    i=0
    while i < imagePixels:
        a = x[i, :]
        aa = np.where(a.mask.astype(int)==0)
        q = np.asarray(aa)
        if q.size < 2:
            i = i+1
        else:
            val.append(np.array([[i, q[0,0]], [i,q[0,q.size-1]]]))
            i = i+1
    cut = np.concatenate(val)
    return cut


## Funtions from Runa ##
def draw_mask(img, cx, cy, r, out=0):
    """ Function to create circular mask
    Created by Runa

    Parameters
    ----------
    img: numpy array
        image to mask
    cx: int [pixel]
        center x of the mask
    cy: int [pixel]
        center y of the mask
    r: int [pixel]
        radius of the mask

    Returns
    -------
    img1: numpy array
        start image mask whit circular new mask
    """
    ss = np.shape(img)
    x = np.arange(ss[0])
    x = np.transpose(np.tile(x, [ss[1], 1]))
    y = np.arange(ss[1])
    y = np.tile(y, [ss[0], 1])
    x = x - cx
    y = y - cy
    nr = np.size(r)
    if nr == 2:
        rr = x*x/r[0]**2+y*y/r[1]**2
        r1 = 1
    else:
        rr = x*x+y*y
        r1 = r**2
    pp = np.where(rr < r1)
    img1 = img.copy()
    if out == 1:
        img1[pp] = 0
    else:
        img1[pp] = 1
    #plt.imshow(img1)
    return img1

def qpupil(mask, xx=None, yy=None, nocircle=0):
    '''
    Function for....
    created by Runa

    Parameters
    ----------
    mask: numpy array

    Returns
    ------
    x0:
    y0:
    r:
    xx: numpy array
        grid of coordinates of the same size as input mask
    yy: numpy array
        grid of coordinates of the same size as input mask
    '''
    idx = np.where(mask == 1)
    ss = np.shape(mask)
    x = np.arange(ss[0]).astype(float)
    x = np.transpose(np.tile(x, [ss[1], 1]))
    y = np.arange(ss[1]).astype(float)
    y = np.tile(y, [ss[0], 1])
    xx = x
    yy = y
    x0 = 0
    y0 = 0
    r = 0
    if nocircle == 0:
        maxv = max(xx[idx])
        minv = min(xx[idx])
        r1 = (maxv-minv)/2
        x0 = r1+minv
        xx = xx - (minv + maxv)/2
        xx = xx/((maxv - minv)/2)
        mx = [minv, maxv]
        maxv = max(yy[idx])
        minv = min(yy[idx])
        r2 = (maxv-minv)/2
        y0 = r2 + minv
        yy = yy - (minv+maxv)/2
        yy = yy/((maxv-minv)/2)
        r = np.mean([r1, r2])
        my = [minv, maxv]
        #imgout  = image[mx[0]:mx[1],my[0]:my[1]]
    return x0, y0, r, xx, yy

def rotate(img, angle):
    ''' Function to rotate the image
    Created by Runa

    Parameters
    ----------
    image: numpy array
        The input array
    angle: float
        The rotation angle in degrees

    Returns
    ------
    img1: numpy array
        The rotated input
    '''
    img1 = ndimage.rotate(img, angle)
    s0 = np.shape(img)
    s1 = np.shape(img1)
    img1 = img1[int((s1[0]-s0[0])/2):s0[0]+int((s1[0]-s0[0])/2),
                int((s1[1]-s0[1])/2):s0[1]+int((s1[1]-s0[1])/2)]
    return img1

def integrate_psd(x,y):
    '''
    to be checked
    '''
    w = np.sqrt(np.sum(y/(2*np.pi*x)**2))
    return w

def crop_frame(mask):
    cir = qpupil(mask)
    cir = np.array(cir[0:3]).astype(int)
    img = mask[cir[0]-cir[2]:cir[0]+cir[2],cir[1]-cir[2]:cir[1]+cir[2]]
    return img


def congrid2D(img, newdims, method='linear', centre=False, minusone=False):
    dims = img.shape
    xo = np.linspace(1,dims[0], dims[0])
    yo = np.linspace(1,dims[1], dims[1])
    x = np.linspace(1,dims[0],newdims[0])
    y = np.linspace(1,dims[1],newdims[1])
    xg,yg =np.meshgrid(x,y,indexing='ij')
    interp = scipy.interpolate.RegularGridInterpolator((xo,yo),img)
    newimg = interp((xg,yg))
    return newimg

def spiral_pos(nstep, step):
    p = np.array([0,0])
    pp = []
    for i in range(nstep):
            direct = (-1)**i*step
            mov = i+1
            for j in range(mov):
                    p = p+np.array([direct,0])
                    pp.append(p)
            for j in range(mov):
                    p = p+np.array([0,direct])
                    pp.append(p)
    pp = np.array(pp)
    plt.plot(pp[:,0],pp[:,1],'-x')
    return pp
