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
from skimage.draw import disk as draw_circle
from m4.ground import zernike
from scipy.ndimage.interpolation import shift
from m4.ground import read_data
from m4.configuration import config_folder_names as config
from matplotlib import pyplot as plt

def patches_analysis(image, radius_m, pixelscale=None, step=None, n_patches=None):
    '''
    Parameters
    ----------
        image: masked array
            image for the analysis
        radius_m: int
            radius of circular patch in meters

    Other Parameters
    ----------
        pixelscale: int
            value of image's pixel scale [px/m]
        step: int
            distance between patches
        n_patches: int
            number of patches for the second cut
            (if it is None sw creates a single crop in the center of the image)
    Returns
    -------
        req: float
            roc at threshold 0.05 or rms at threshold 0.95
        list_ima: list
            list of images used for the analysis
        result_vect: numpy array
            vector containing the test analysis results
    '''
    if pixelscale is None:
        ps = 1/OttParameters.pscale
    else:
        ps = 1/pixelscale

    print(ps)
    raggio_px = radius_m/ps

    nn = image.compressed().shape[0]
    idx = np.where(image.mask==0)
    x = idx[0]
    y = idx[1]

    if step is None:
        n_point = nn
        step = 1
    else:
        n_point = np.int(nn/step)+1

    result_list = []
    list_ima = []
    list_big = []
    for i in range(n_point):
        p = i + i * (step - 1)
        #print('%d' %p)
        #print('%d %d' %(x[p], y[p]))
        if radius_m == 0.04:
            thresh = 0.05
            ima = _circleImage(image, x[p], y[p], raggio_px)
            if ima is not None:
                #list_ima.append(ima)
                alpha, beta = curv_fit_v2(ima, 1/ps)
                raggi_km = roc(alpha, beta)
                #raggio = roc(2*radius_m, alpha, beta)
                #print(raggio)
                result_list.append(raggi_km)
        if radius_m == 0.015:
            thresh = 0.95
            r1 = 0.1/2
            r_px1 = r1/ps
            ima = _circleImage(image, x[p], y[p], r_px1)
            if ima is not None:
                new_ima = tiptilt_fit(ima)
                #list_big.append(new_ima)
                if n_patches is None:
                    final_ima = _circleImage(new_ima, x[p], y[p], raggio_px)
                    if final_ima is not None:
                           if i == 2:
                            list_ima.append(final_ima)
                        final_ima = final_ima - np.mean(final_ima)
                        result_list.append(np.std(final_ima))
                else:
                    list_circleima = _circleImageList(new_ima, x[p], y[p], n_patches, raggio_px)
                    if list_circleima is not None:
                        for ima in list_circleima:
                            list_ima.append(ima)
                            result_list.append(np.std(ima))
        if radius_m == 0.25:
            thresh = 0.95
            r_px1 = 0.1/2/ps
            ima = _circleImage(image, x[p], y[p], r_px1)
            if ima is not None:
                new_ima = tiptilt_fit(ima)
                #list_ima.append(new_ima)
                result_list.append(np.std(new_ima))

    result_vect = np.array(result_list)
    result_sort = np.copy(result_vect)
    result_sort.sort()
    dim = result_sort.size
    req = result_sort[np.int(thresh*dim)]
    return req, list_ima, result_vect


def _circleImageList(new_ima, x, y, n_point, r_px):
    '''
    Parameters
    ----------
        image: masked array
            image for the analysis
        x: int
            coordinate x
        y: int
            coordinate y
        n_point: int
            number of images for second cut
        r_px: int
            radius of circular patch in meters
    Returns
    -------
        final_ima_list: list
            list of circle cropped image
    '''
    if new_ima is not None:
        nn = new_ima.compressed().shape[0]
        idx = np.where(new_ima.mask==0)
        x = idx[0]
        y = idx[1]
        step = np.int(nn/(n_point-1))
        final_ima_list = []
        for i in range(n_point):
            p = i + i * (step - 2)
            final_ima = _circleImage(new_ima, x[p], y[p], r_px)
            if final_ima is not None:
                final_ima_list.append(final_ima)
        return final_ima_list

def _circleImage(image, x, y, raggio_px):
    ''' Function to create circular cuts of the image.
    The function excludes cuts that come out of the edges of the image
    and check the threshold for valid point
    NOTE: if image is out of conditions return None type

    Parameters
    ----------
        image: masked array
            image for the analysis
        x: int
            coordinate x
        y: int
            coordinate y
        radius_px: int
            radius of circular patch in pixels

    Returns
    -------
        ima: numpy masked array
            cut image
    '''
    th_validPoint = np.pi * raggio_px**2 * 30 / 100
    if image is not None:
        circle = np.ones((image.shape[0],image.shape[1])).astype(int)
        rr, cc = draw_circle((x, y), raggio_px)
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
            valid_point = ima.compressed().shape[0]
            if valid_point >= th_validPoint:
                return ima

def tiptilt_fit(ima):
    '''
    Parameters
    ----------
        image: masked array

    Returns
    -------
        ima_ttr: numpy masked array
            image without tip and tilt
    '''
    if ima is not None:
        coef, mat = zernike.zernikeFit(ima,
                                    np.array([2, 3]))
        surf = zernike.zernikeSurface(ima, coef, mat)
        new_image = ima - surf
        ima_ttr = np.ma.masked_array(new_image, mask=ima.mask)
        return ima_ttr

# def curv_fit(image, test_diameter):
#     '''
#     Parameters
#     ----------
#         image: masked array
#             image for the analysis
#         test_diameter: int
#             diameter for xy coordinates
#
#     Returns
#     -------
#         alpha: float
#             analytical coefficient of scalloping
#         beta: float
#             analytical coefficient of scalloping
#     '''
#     size = np.array([image.shape[0], image.shape[1]])
#     ima_x = np.arange(size[0], dtype = float)
#     ima_y = np.arange(size[1], dtype = float)
#     xx = np.tile(ima_x, (size[0], 1))
#     yy = np.tile(ima_y, (size[1], 1)).T
#
#     idx = np.where(image.mask==0)
# #    xx = xx*(test_diameter/2)
# #    yy = yy*(test_diameter/2)
#
#     xx = xx*(test_diameter)/size[0]
#     yy = yy*(test_diameter)/size[1]
#
#     nn = image.compressed().shape[0]
#     zmat = np.zeros((nn, 6))
#     zmat[:,0] = xx[idx]**2
#     zmat[:,1] = yy[idx]**2
#     zmat[:,2] = xx[idx]*yy[idx]
#     zmat[:,3] = xx[idx]
#     zmat[:,4] = yy[idx]
#     zmat[:,5] = np.ones(nn)
#
#     inv = np.linalg.pinv(zmat)
#     coeff = np.dot(inv, image.compressed())
#
#     alpha = 0.5*( coeff[0]+ coeff[1]+ np.sqrt((coeff[0]-coeff[1])**2 + coeff[2]**2) )
#     beta  = 0.5*( coeff[0]+ coeff[1]- np.sqrt((coeff[0]-coeff[1])**2 + coeff[2]**2) )
#     return alpha, beta

def curv_fit_v2(image, platescale_px_mm):
    '''
    Parameters
    ----------
        image: masked array
            image for the analysis
        platescale_px_mm: double
            platescale in pixel / mm

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
    xx = xx/platescale_px_mm
    yy = yy/platescale_px_mm

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

    #compute the curvatures (1/R) and return them in km
    alpha = 0.5*( coeff[0]+ coeff[1]+ np.sqrt((coeff[0]-coeff[1])**2 + coeff[2]**2) )
    beta  = 0.5*( coeff[0]+ coeff[1]- np.sqrt((coeff[0]-coeff[1])**2 + coeff[2]**2) )
    return alpha, beta

def roc(alpha, beta):
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
    vect = np.array([1/(2*alpha), 1/(2*beta)])
    raggi_km = np.abs(vect)/1000
#     wfe = test_diameter**2 /(8*np.sqrt(3)) * np.sqrt(2*(alpha-beta)**2 - (alpha+beta)**2)
#     rho = test_diameter/2
#     raggio = (rho**2 + wfe**2)/(2*wfe)
    return raggi_km

def slope(image, pscale):
    '''
    Parameters
    ----------
        image: masked array
            image for the analysis
        pscale: float
            pixel scale of image [px/m]
    Returns
    -------
        slope: numpy masked array
    '''
    ax = ((image - shift(image, [1,0]))*pscale)
    ay = ((image - shift(image, [0,1]))*pscale)
    mask = np.ma.mask_or(image.mask, shift(image.mask, [1,0]))
    sp = np.sqrt(((ax))**2+((ay))**2)
    #s = np.sqrt((np.arctan(ax))**2+(np.arctan(ay))**2)
    masked_slope = np.ma.masked_array(sp, mask=mask)
    return masked_slope


### REQ ###

def test242(image, pscale=None):
    '''
    Parameters
    ----------
        image: masked array
            robust image for the analysis

    Returns
    -------
        rms: float
            rms slope in arcsec
    '''
    mask = np.invert(image.mask).astype(int)
    x, y, r, xx, yy = geo.qpupil(mask)
    if pscale is None:
        pscale = r * (1/0.17)
    else:
        pscale = pscale

    print(pscale)
    sp = slope(image, pscale) #pscale [px/m]
    # sp in pixel * fattore di conversione da rad ad arcsec
    slope_arcsec = sp * 206265
    rms = slope_arcsec.std()
    return rms

def test243(image, radius_m, pscale=None, step=None, n_patches=None):
    ''' Return rms at the interactuator scale 31 mm or 150 mm (thresh = 0.95)

    Parameters
    ----------
        image: masked array
            robust image for the analysis
        radius_m: int
            radius of circular patch in meters

    Other Parameters
    ----------------
        step: int
            distance between patches
        n_patches: int
            number of patches for the second cut
            (if it is None sw creates a single crop in the center of the image)

    Returns
    -------
        rms: float
            rms at threshold 0.95
'''
    # radius_m = 0.015 or 0.1
    mask = np.invert(image.mask).astype(int)
    x, y, r, xx, yy = geo.qpupil(mask)
    if pscale is None:
        pscale = r * (1/0.17)
    else:
        pscale = pscale

    rms, list_ima, result_vect = patches_analysis(image, radius_m, pscale, step, n_patches)
    return rms

def test283(image, pscale=None, step=None):
    ''' Return roc on 80 mm spatial scale (thresh = 0.05)

    Parameters
    ----------
        image: masked array
            robust image for the analysis

    Other Parameters
    ----------------
        step: int
            distance between patches

    Returns
    -------
        roc: float
            roc at threshold 0.05
'''
    mask = np.invert(image.mask).astype(int)
    x, y, r, xx, yy = geo.qpupil(mask)
    if pscale is None:
        pscale = r * (1/0.17) #pscale[px/m]
    else:
        pscale = pscale

    print(pscale)
    roc, list_ima, result_vect = patches_analysis(image, 0.04, pscale, step)
    return roc

def diffPiston(image):
    '''
    Parameters
    ----------
        image: masked array
            image for the analysis

    Returns
    -------
        diff_piston: numpy masked array
    '''
    dimx = image.shape[0]
    dimy = image.shape[1]
    imas = image[:, 0:np.int(dimy/2)]
    imad = image[np.int(dimx/2):, :]

    coefs, mat = zernike.zernikeFit(imas, np.arange(3)+1)
    coefd, mat = zernike.zernikeFit(imad, np.arange(3)+1)
    diff_piston = coefs[0]-coefd[1]
    return diff_piston


### ROBUST IMAGE ###
def imageOpticOffset(data_file_path, start, stop):
    '''
    Parameters
    ----------
    data_file_path: string
        data file path for measurement to analyze
    start: int
        number of first image to use for the data analysis
    stop: int
        last number of measurement to use

    Returns
    -------
    image: numpy  masked array
        mean image of the selected data
    '''
    last_name = data_file_path.split('/')[-1]
    if last_name == 'hdf5':
        list = glob.glob(os.path.join(data_file_path, '*.h5'))
        tt = data_file_path.split('/')[-2]
        #ext = 1
    else:
        list = glob.glob(os.path.join(data_file_path, '*.fits'))
        tt = data_file_path.split('/')[-1]
        #ext = 0

    list.sort()
    list = list[start:stop]

    cube = None
    print('Creating cube for offset image:')
    for name in list:
        #nn = name.split('/')[-1]
        #print(nn)
        image = read_data.read_phasemap(name)
        if cube is None:
            cube = image
        else:
            cube = np.ma.dstack((cube, image))

    image = np.mean(cube, axis=2)

    results_path = os.path.join(config.OUT_FOLDER, 'Req', tt)
    fits_file_name = os.path.join(results_path, 'OptOffset.fits')
    pyfits.writeto(fits_file_name, image.data, overwrite=True)
    pyfits.append(fits_file_name, image.mask.astype(int), overwrite=True)
    return image


def robustImageFromDataSet(n_images, data_file_path, zernike_vector_to_subtract, offset=None):
    ''' From fits files and whit offset subtraction

    Parameters
    ----------
    n_images: int
        number of images to analyze
    path: string
        total path for data analysis

    Other Parameters
    ----------------
    offset: if it is None data analysis is made by split n_images in two
            else re-reads the offset image saved in the tt folder and subtracts it
            to each image during cube creation

    Returns
    -------
    robust_image: numpy masked array
        robust image from data set
    '''
    last_name = data_file_path.split('/')[-1]
    if last_name == 'hdf5':
        list_tot = glob.glob(os.path.join(data_file_path, '*.h5'))
        tt = data_file_path.split('/')[-2]
        #ext = 1
    else:
        list_tot = glob.glob(os.path.join(data_file_path, '*.fits'))
        tt = data_file_path.split('/')[-1]
        #ext = 0

    list_tot.sort()
    list = list_tot[0: n_images]
    if offset is None:
        half = np.int(len(list)/2)
        list1 = list[0:half]
        list2 = list[half:]

        cube1 = None
        print('Creating cube 1')
        for name in list1:
            #print(name)
            image = read_data.read_phasemap(name)
            if cube1 is None:
                cube1 = image
            else:
                cube1 = np.ma.dstack((cube1, image))

        cube2 = None
        print('Creating cube 2')
        for name in list2:
            #print(name)
            image = read_data.read_phasemap(name)
            if cube2 is None:
                cube2 = image
            else:
                cube2 = np.ma.dstack((cube1, image))

        mean1 = np.ma.mean(cube1, axis=2)
        mean2 = np.ma.mean(cube2, axis=2)

        final_image = mean2 - mean1

    else:
        fits_file_name = os.path.join(config.OUT_FOLDER, 'Req', tt,
                                      'OptOffset.fits')
        image_optOffset = read_data.readFits_maskedImage(fits_file_name)

        cube = None
        print('Creating cube')
        for name in list:
            #print(name)
            ima = read_data.read_phasemap(name)
            image = ima - image_optOffset
            if cube is None:
                cube = image
            else:
                cube = np.ma.dstack((cube, image))
        final_image = np.ma.mean(cube, axis=2)

    coef, mat = zernike.zernikeFit(final_image, zernike_vector_to_subtract)
    surf = zernike.zernikeSurface(final_image, coef, mat)
    image_ttr = final_image - surf
    return image_ttr





###TEST###
def _imaTest():
    ff = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/ZernikeCommandTest/20191210_110019/mode0002_measure_segment00_neg.fits'
    hduList=pyfits.open(ff)
    image = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
    return image

#    path_rm = '/Users/rm/Desktop/Arcetri/M4/RM_20201020.fits'
#    path_par = '/Users/rm/Desktop/Arcetri/M4/PM_20191113.fits'
def _readTestData(path):
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
