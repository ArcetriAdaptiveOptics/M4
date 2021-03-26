'''
Authors
  - C. Selmi:  written in October 2020
'''

import numpy as np
import os
import glob
from m4.ground import geo
from m4.configuration.ott_parameters import OttParameters
from astropy.io import fits as pyfits
from skimage.draw import circle as draw_circle
from m4.ground import zernike
from scipy.ndimage.interpolation import shift
from m4.ground import read_data
from m4.configuration import config

def patches_analysis(image, radius_m, fit, pixelscale=None, step=None, n_patches=None):
    '''
    Parameters
    ----------
        image: masked array
            image for the analysis
        radius_m: int
            radius of circular patch in meters
        fit: int
            0 to remove tip e tilt
            other values to not remove it
    Other Parameters
    ----------
        pixelscale: int
            value of image's pixel scale
        step: int
            distance between patches
        n_patches: int
        number of patches for the second cut
        (if it is None sw creates a single crop in the center of the image)
    Returns
    -------
        rms: numpy array
            std of every patches
        ord_rms: numpy array
            tidy std of every patches
    '''
    diameter = image.shape[0]
    #OttParameters.PARABOLA_PUPIL_XYRADIUS[2] = np.int(diameter/2)
    #metri = 0.05
    if pixelscale is None:
        ps = 1/OttParameters.pscale
    else:
        ps = pixelscale

    raggio_px = radius_m/ps

    nn = image.compressed().shape[0]
    th_validPoint = np.pi * raggio_px**2 * 30 / 100
    idx = np.where(image.mask==0)
    x = idx[0]
    y = idx[1]

    if step is None:
        n_point = nn
        step = 1
    else:
        n_point = np.int(nn/step)+1

    rms_list = []
    list_ima = []
    for i in range(n_point):
        p = i + i * (step - 1)
        print('%d' %p)
#         circle = np.ones((image.shape[0],image.shape[1])).astype(int)
#         rr, cc = draw_circle(x[p], y[p], raggio_px)
#         circle[rr, cc] = 0
#         mask_prod = np.ma.mask_or(circle.astype(bool), image.mask)
#         ima = np.ma.masked_array(image, mask=mask_prod)
        print('%d %d' %(x[p], y[p]))
        if n_patches is None:
            ima = _patchesAndFit(image, radius_m, x[p], y[p], fit, ps)
            if ima is None:
                pass
            else:
                list_ima.append(ima)
                valid_point = ima.compressed().shape[0]
                if valid_point >= th_validPoint:
                    rms_list.append(np.std(ima))
                else:
                    pass
        else:
            list_ima = _patchesAndFitMatrix(image, radius_m, x[p], y[p], n_patches, ps)
            if list_ima is None:
                pass
            else:
                for ima in list_ima:
                    valid_point = ima.compressed().shape[0]
                    if valid_point >= th_validPoint:
                        rms_list.append(np.std(ima))
                    else:
                        pass
    rms = np.array(rms_list)
    ord_rms = np.sort(rms)
    return rms, ord_rms

def _patchesAndFit(image, radius_m, x, y, fit, pixelscale=None):
    '''
    Parameters
    ----------
        image: masked array
            image for the analysis
        radius_m: int
            radius of circular patch in meters
        x: int
            coordinate x
        y: int
            coordinate y
        fit: int
            0 to remove tip e tilt
            other values to not remove it
    Other Parameters
    ----------
        pixelscale: int
            value of image's pixel scale
    Returns
    -------
        final_ima: numpy masked array
            circle cropped image
    '''
    if pixelscale is None:
        ps = 1/OttParameters.pscale
    else:
        ps = pixelscale

    if fit == 0:
        if radius_m<=(0.1/2):
            r1 = 0.1/2
            r_px1 = r1/ps
            ima = _circleImage(image, x, y, r_px1)
            new_ima = surf_fit(ima)
            r2 = radius_m
            r_px2 = r2/ps
            final_ima = _circleImage(new_ima, x, y, r_px2)
        else:
            r_px2 = radius_m/ps
            ima = _circleImage(image, x, y, r_px2)
            final_ima = surf_fit(ima)
    else:
        r1 = radius_m
        r_px1 = r1/ps
        final_ima = _circleImage(image, x, y, r_px1)

    return final_ima

def _patchesAndFitMatrix(image, radius_m, x, y, n_point, pixelscale=None):
    '''
    Parameters
    ----------
        image: masked array
            image for the analysis
        radius_m: int
            radius of circular patch in meters
        x: int
            coordinate x
        y: int
            coordinate y
        n_point: int
            number of images for second cut
    Other Parameters
    ----------
        pixelscale: int
            value of image's pixel scale
    Returns
    -------
        final_ima_list: list
            list of circle cropped image
    '''
    if pixelscale is None:
        ps = 1/OttParameters.pscale
    else:
        ps = pixelscale

    r1 = 0.1/2
    r_px1 = r1/ps
    ima = _circleImage(image, x, y, r_px1)
    new_ima = surf_fit(ima)
    r2 = radius_m
    r_px2 = r2/ps
    if new_ima is None:
        pass
    else:
        nn = new_ima.compressed().shape[0]
        idx = np.where(new_ima.mask==0)
        x = idx[0]
        y = idx[1]
        step = np.int(nn/(n_point-1))
        final_ima_list = []
        for i in range(n_point):
            p = i + i * (step - 2)
            final_ima = _circleImage(new_ima, x[p], y[p], r_px2)
            if final_ima is None:
                pass
            else:
                final_ima_list.append(final_ima)
        return final_ima_list

def _circleImage(image, x, y, raggio_px):
    ''' Function to create circular cuts of the image.
    The function excludes cuts that come out of the edges of the image
    '''
    if image is None:
        pass
    else:
        circle = np.ones((image.shape[0],image.shape[1])).astype(int)
        rr, cc = draw_circle(x, y, raggio_px)
        test_r = (len(list(filter (lambda x : x >= image.shape[0], rr))) > 0)
        test_r2 = (len(list(filter (lambda x : x <= 0, rr))) > 0)
        test_c = (len(list(filter (lambda x : x >= image.shape[0], cc))) > 0)
        if test_r == True:
            pass
        elif test_r2 == True:
            pass
        elif test_c == True:
            pass
        else:
            circle[rr, cc] = 0
            mask_prod = np.ma.mask_or(circle.astype(bool), image.mask)
            ima = np.ma.masked_array(image, mask=mask_prod)
            return ima

def surf_fit(ima):
    '''
    Parameters
    ----------
        image: masked array

    Returns
    -------
        pp: numpy masked array
            image without tip and tilt
    '''
    if ima is None:
        pass
    else:
        coef, mat = zernike.zernikeFit(ima,
                                    np.array([2, 3]))
        new_image = zernike.zernikeSurface(ima, coef, mat)
        pp = np.ma.masked_array(new_image, mask=ima.mask) 
        return pp

def curv_fit(image, test_diameter):
    '''
    Parameters
    ----------
        image: masked array
            image for the analysis
        test_diameter: int
            diameter for xy coordinates

    Returns
    -------
        alpha: float
            analytical coefficient of scalloping
        beta: float
            analytical coefficient of scalloping
    '''
    size = np.array([image.shape[0], image.shape[1]])
    ima_x = np.arange(size[0], dtype = float)
    ima_y = np.arange(size[1], dtype = float)
    xx = np.tile(ima_x, (size[0], 1))
    yy = np.tile(ima_y, (size[1], 1)).T

    idx = np.where(image.mask==0)
    xx = xx*(test_diameter/2)
    yy = yy*(test_diameter/2)

    nn = image.compressed().shape[0]
    zmat = np.zeros((nn, 6))
    zmat[:,0] = xx[idx]**2
    zmat[:,1] = yy[idx]**2
    zmat[:,2] = xx[idx]*yy[idx]
    zmat[:,3] = xx[idx]
    zmat[:,4] = yy[idx]
    zmat[:,5] = np.ones(nn)

    inv = np.linalg.pinv(zmat)
    coeff = np.dot(inv, image.compressed())

    alpha = 0.5*( coeff[1]+ coeff[2]+ np.sqrt((coeff[1]-coeff[2])**2 + coeff[3]**2) )
    beta  = 0.5*( coeff[1]+ coeff[2]- np.sqrt((coeff[1]-coeff[2])**2 + coeff[3]**2) )
    return alpha, beta

def roc(test_diameter, alpha, beta):
    '''
    Parameters
    ----------
        test_diameter: int
            diameter for xy coordinates
        alpha: float
            analytical coefficient of scalloping
        beta: float
            analytical coefficient of scalloping

    Returns
    -------
        raggio: float
            radius of curvature
    '''
    wfe = test_diameter**2 /(8*np.sqrt(3)) * np.sqrt(2*(alpha-beta)**2 - (alpha+beta)**2)
    rho = test_diameter/2
    raggio = (rho**2 + wfe**2)/(2*wfe)
    return raggio

def slope(image):
    '''
    Parameters
    ----------
        image: masked array
            image for the analysis
    Returns
    -------
        slope: masked array
    '''
    mask = np.invert(image.mask).astype(int)
    x, y, r, xx, yy = geo.qpupil(mask)
    pscale = r * (1/0.17)

    ax = ((image - shift(image, [1,0]))*pscale)
    ay = ((image - shift(image, [0,1]))*pscale)
    mask = np.ma.mask_or(image.mask, shift(image.mask, [1,0]))
    sp = np.sqrt(((ax))**2+((ay))**2)
    #s = np.sqrt((np.arctan(ax))**2+(np.arctan(ay))**2)
    masked_slope = np.ma.masked_array(sp, mask=mask)
    return masked_slope


### REQ ###

def test242(image):
    sp = slope(image)
    # sp in pixel * fattore di conversione da rad ad arcsec
    slope_arcsec = sp * 206265
    rms = slope_arcsec.std()
    return rms

def test243(image, ps):
    diameter = 0 #average inter-actuator spacing
    fit = 0
    step = None
    rms, rms_ord = patches_analysis(image, diameter/2, fit, ps, step)
    return rms

def test283(image):
    diameter = 0.08
    alpha, beta = curv_fit(image, diameter)
    raggio = roc(diameter, alpha, beta)
    return raggio

def diffPiston(image):
    dimx = image.shape[0]
    dimy = image.shape[1]
    imas = image[:, 0:np.int(dimy/2)]
    imad = image[np.int(dimx/2):, :]

    coefs, mat = zernike.zernikeFit(imas, np.arange(3)+1)
    coefd, mat = zernike.zernikeFit(imad, np.arange(3)+1)
    diff_piston = coefs[0]-coefd[1]
    return diff_piston


### ROBUST IMAGE ###

def imageOpticOffset(data_file_path):
    list = glob.glob(os.path.join(data_file_path, '*.fits'))
    list.sort()

    cube = None
    print('Creating cube for offset image:')
    for name in list:
        nn = name.split('/')[-1]
        print(nn)
        image = read_data.readFits_maskedImage(name)
        if cube is None:
            cube = image
        else:
            cube = np.ma.dstack((cube, image))

    image = np.mean(cube, axis=2)
    dove = '/home/labot/data/M4/Data/M4Data/Results'
    fits_file_name = os.path.join(dove, 'OptOffset.fits')
    pyfits.writeto(fits_file_name, image.data)
    pyfits.append(fits_file_name, image.mask.astype(int))
    return image

def robustImageFromDataset(path):
    ''' From h5 files '''
    list = glob.glob(os.path.join(path, '*.h5'))
    list.sort()
    half = np.int(len(list)/2)
    list1 = list[0:half]
    list2 = list[half:]

    cube1 = None
    print('Creating cube 1:')
    for name in list1:
        print(name)
        image = read_data.InterferometerConverter.from4D(name)
        if cube1 is None:
            cube1 = image
        else:
            cube1 = np.ma.dstack((cube1, image))

    cube2 = None
    print('Creating cube 2:')
    for name in list2:
        print(name)
        image = read_data.InterferometerConverter.from4D(name)
        if cube2 is None:
            cube2 = image
        else:
            cube2 = np.ma.dstack((cube1, image))

    mean1 = np.ma.mean(cube1, axis=2)
    mean2 = np.ma.mean(cube2, axis=2)

    image = mean2 -mean1
    return image

def robustImageFromStabilityData(n_images, path):
    ''' From fits files and whit offset subtraction'''
    list_tot = glob.glob(os.path.join(path, '*.fits'))
    list_tot.sort()
    list = list_tot[0:n_images]
    half = np.int(len(list)/2)
    list1 = list[0:half]
    list2 = list[half:]

    cube1 = None
    print('Creating cube 1:')
    for name in list1:
        print(name)
        image = read_data.readFits_maskedImage(name)
        if cube1 is None:
            cube1 = image
        else:
            cube1 = np.ma.dstack((cube1, image))

    cube2 = None
    print('Creating cube 2:')
    for name in list2:
        print(name)
        image = read_data.readFits_maskedImage(name)
        if cube2 is None:
            cube2 = image
        else:
            cube2 = np.ma.dstack((cube1, image))

    mean1 = np.ma.mean(cube1, axis=2)
    mean2 = np.ma.mean(cube2, axis=2)

    image = mean2 -mean1

    fits_file_name = os.path.join(config.path_name.OUT_FOLDER, 'Req', 'OptOffset5000.fits')
    image_optOffset = read_data.readFits_maskedImage(fits_file_name)
    final_image = image - image_optOffset
    return final_image





###TEST###
def imaTest():
    ff = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/ZernikeCommandTest/20191210_110019/mode0002_measure_segment00_neg.fits'
    hduList=pyfits.open(ff)
    image = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
    return image

#    path_rm = '/Users/rm/Desktop/Arcetri/M4/RM_20201020.fits'
#    path_par = '/Users/rm/Desktop/Arcetri/M4/PM_20191113.fits'
def readTestData(path):
    hduList = pyfits.open(path)
    ogg = hduList[0].data
    dim = int(np.sqrt(ogg[:,0].size))
#     image = np.zeros((dim, dim))
#     xx = np.reshape(ogg[:,1], [dim,dim])
#     yy = np.reshape(ogg[:,0], [dim,dim])
#     z = ogg[:,2]
    zz = np.reshape(ogg[:,2], [dim,dim])
    mask = np.isnan(zz)
    prova = np.ma.masked_array(zz, mask=mask)
    ps = ogg[:,1][1] - ogg[:,1][0]
    return prova, ps

### Per l'immagine 591X591
# vv = np.ma.masked_array(np.zeros((1,image.shape[0])), mask=np.ones((1,image.shape[0])).astype(bool))
# vv2 =np.ma.masked_array(np.zeros((image.shape[0]+1,1)), mask=np.ones((image.shape[0]+1, 1)).astype(bool)) 
# pp = np.ma.append(image, vv, axis=0)
# new = np.ma.append(pp, vv2, axis=1)
