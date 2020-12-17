# -*- coding: utf-8 -*-
"""
Autors
  - R. Briguglio: created on Mon Mar 16 11:00:08 2020

"""
import numpy as np
#import image as image
#import scipy
from matplotlib import pyplot as plt
from scipy import ndimage

def draw_mask(img, cx, cy, r, out=0):
    """ Function to create circular mask
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
        y0 = r2 +minv
        yy = yy - (minv+maxv)/2
        yy = yy/((maxv-minv)/2)
        r = np.mean([r1, r2])
        my = [minv, maxv]
        #imgout  = image[mx[0]:mx[1],my[0]:my[1]]
    return x0, y0, r, xx, yy

def rotate(img, angle):
    ''' Function to rotate the image
    '''
    img1 = ndimage.rotate(img, angle)
    s0 = np.shape(img)
    s1 = np.shape(img1)
    img1 = img1[int((s1[0]-s0[0])/2):s0[0]+int((s1[0]-s0[0])/2),
                int((s1[1]-s0[1])/2):s0[1]+int((s1[1]-s0[1])/2)]
    return img1
