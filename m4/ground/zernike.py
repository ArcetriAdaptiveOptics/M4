#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

"""
@file zernike.py
@brief Zernike generation library
@author M.Xompero
@url -
@date 20200202

Created by Tim van Werkhoven (t.i.m.vanwerkhoven@xs4all.nl) on 2011-10-12
Copyright (c) 2011 Tim van Werkhoven. All rights reserved.

This file is licensed under the Creative Commons Attribution-Share Alike
license versions 3.0 or higher, see
http://creativecommons.org/licenses/by-sa/3.0/

HOW TO USE IT::

    from m4.ground import zernike
    coeff, mat = zernike.zernikeFit(img, zernike_index_vector)
    surf_image = zernike.zernikeSurface(img, coef, mat)
"""

### Libraries

import numpy as np
from m4.ground import geo
fac = np.math.factorial

def zernikeFit(img, zernike_index_vector):
    '''
    Parameters
    ----------
    img: numpy masked array
        image for zernike fit
    zernike_index_vector: numpy array
        vector containing the index of Zernike modes to be fitted starting from 1

    Returns
    -------
    coeff: numpy array [m]
        vector of zernike coefficients
    mat: numpy array
    '''
    img1 = img.data
    mask = np.invert(img.mask).astype(int)
    x, y, r, xx, yy = geo.qpupil_circle(img)
    mm = (mask==1)
    coeff = _surf_fit(xx[mm], yy[mm], img1[mm], zernike_index_vector)
    mat = _getZernike(xx[mm], yy[mm], zernike_index_vector)
    return coeff, mat

def zernikeFitAuxmask(img, auxmask,zernike_index_vector):
    '''
    Parameters
    ----------
    img: numpy masked array
        image for zernike fit
    zernike_index_vector: numpy array
        vector containing the index of Zernike modes to be fitted starting from 1

    Returns
    -------
    coeff: numpy array [m]
        vector of zernike coefficients
    mat: numpy array
    '''
    img1 = img.data
    mask = np.invert(img.mask).astype(int)
    x, y, r, xx, yy = geo.qpupil(auxmask)
    mm = (mask==1)
    coeff = _surf_fit(xx[mm], yy[mm], img1[mm], zernike_index_vector)
    mat = _getZernike(xx[mm], yy[mm], zernike_index_vector)
    return coeff, mat



def zernikeSurface(img, coef, mat):
    '''
    Parameters
    ----------
    img: numpy masked array
        image for zernike fit
    coeff: numpy array [m]
        vector of zernike coefficients
    mat: numpy array

    Returns
    -------
    surf: numpy masked array
        zernike surface generate by coeff
    '''
#     img1 = img.data
#     mask = np.invert(img.mask).astype(int)
#     x, y, r, xx, yy = geo.qpupil(mask)
#     mm = (mask==1)
#     aa = getZernike(xx[mm], yy[mm], zernike_index_vector)
#     coeff = surf_fit(xx[mm], yy[mm], img1[mm], zernike_index_vector)
    mm = np.where(img.mask == 0)
    zernike_surface = np.zeros((img.shape[0], img.shape[1]))
    zernike_surface[mm] = np.dot(mat, coef)
    surf = np.ma.masked_array(zernike_surface, mask=img.mask)
    return surf

def _surf_fit(xx, yy, zz, zlist, ordering='noll'):
    A = _getZernike(xx, yy, zlist, ordering)
    B = np.transpose(zz.copy())
    coeff = (np.linalg.lstsq(A, B, rcond=-1))[0]
    return coeff


### Init functions
def _getZernike(xx,yy,zlist,ordering='noll'):
    if min(zlist) ==0:
        #print("Zernike index must be greater or equal to 1")
        raise OSError("Zernike index must be greater or equal to 1")
    rho = np.sqrt(yy**2 + xx**2)
    phi = np.arctan2(yy, xx)

    #rho /= scale_factor
    zkm = []
    norm = []
    #for k, j in enumerate([0, 1, 2, 4]):
    for j in zlist:

        if ordering=='noll':
            m, n = _l2mn_noll(j)
            cnorm = np.sqrt(n+1) if m == 0 else np.sqrt(2.0*(n+1))
        elif ordering=='ansi': #da rivedere ordine e normalizzazione
            m, n = _l2mn_ansi(j)
            cnorm = 1
        zkm.append(cnorm*_zernike(m, n, rho, phi))
    return np.transpose(np.array(zkm))

def _zernike_rad(m, n, rho):
    """
    Calculate the radial component of Zernike polynomial (m, n)
    given a grid of radial coordinates rho.

    >>> zernike_rad(3, 3, 0.333)
    0.036926037000000009
    >>> zernike_rad(1, 3, 0.333)
    -0.55522188900000002
    >>> zernike_rad(3, 5, 0.12345)
    -0.007382104685237683
    """

    if (n < 0 or m < 0 or abs(m) > n):
        raise ValueError

    if ((n-m) % 2):
        return rho*0.0

    pre_fac = lambda k: (-1.0)**k * fac(n-k) / ( fac(k) * fac( (n+m)/2.0 - k ) * fac( (n-m)/2.0 - k ) )

    return sum(pre_fac(k) * rho**(n-2.0*k) for k in range((n-m)//2+1))

def _zernike(m, n, rho, phi):
    """
    Calculate Zernike polynomial (m, n) given a grid of radial
    coordinates rho and azimuthal coordinates phi.

    >>> zernike(3,5, 0.12345, 1.0)
    0.0073082282475042991
    >>> zernike(1, 3, 0.333, 5.0)
    -0.15749545445076085
    """
    if (m > 0): return _zernike_rad(m, n, rho) * np.cos(m * phi)
    if (m < 0): return _zernike_rad(-m, n, rho) * np.sin(-m * phi)
    return _zernike_rad(0, n, rho)

def _zernikel(j, rho, phi):
    """
    Calculate Zernike polynomial with Noll coordinate j given a grid of radial
    coordinates rho and azimuthal coordinates phi.

    >>> zernikel(0, 0.12345, 0.231)
    1.0
    >>> zernikel(1, 0.12345, 0.231)
    0.028264010304937772
    >>> zernikel(6, 0.12345, 0.231)
    0.0012019069816780774
    """
    n = 0
    while (j > n):
        n += 1
        j -= n

    m = -n+2*j
    return _zernike(m, n, rho, phi)

def _l2mn_ansi(j):
    n = 0
    while (j > n):
        n += 1
        j -= n

    m = -n+2*j

    return m, n

def _l2mn_noll(j):
    """
    Find the [n,m] list giving the radial order n and azimuthal order
    of the Zernike polynomial of Noll index j.

    Parameters:
        j (int): The Noll index for Zernike polynomials

    Returns:
        list: n, m values
    """
    n = int((-1.+np.sqrt(8*(j-1)+1))/2.)
    p = (j-(n*(n+1))/2.)
    k = n%2
    m = int((p+k)/2.)*2 - k

    if m!=0:
        if j%2==0:
            s=1
        else:
            s=-1
        m *= s

    return [m, n]


##### TEST

def _test():
    from matplotlib.pyplot import imshow,show
    import numpy as np
    x = np.arange(-1,1,0.01)
    yy = np.tile(np.transpose(x),(200,1))
    xx = np.transpose(yy)
    mm = (xx**2+yy**2) < 1
    aa = _getZernike(xx[mm],yy[mm],np.arange(11)+1);
    ii=xx*0;ii[mm]=aa[:,9];imshow(ii);show()


