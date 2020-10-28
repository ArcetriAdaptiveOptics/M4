'''
Autors
  - C. Selmi:  written in October 2020
'''

from m4.configuration.ott_parameters import *
from astropy.io import fits as pyfits
from skimage.draw import circle as draw_circle
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from scipy.ndimage.interpolation import shift

def patches_extractor(image, radius_m, step=None):
    '''
    Parameters
    ----------
        image: masked array
            image for the analysis
        radius_m: int
            radius of circular patch in meters
    Other Parameters
    ----------
        step: int
            distance between patches
    Returns
    -------
        rms: numpy array
            std of every patches
        ord_rms: numpy array
            tidy std of every patches
    '''
    #metri = 0.05
    ps = 1/OttParameters.pscale
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

    rms = np.zeros(n_point)
    for i in range(n_point):
        p = i + i * (step - 1)
        print('%d' %p)
#         circle = np.ones((image.shape[0],image.shape[1])).astype(int)
#         rr, cc = draw_circle(x[p], y[p], raggio_px)
#         circle[rr, cc] = 0
#         mask_prod = np.ma.mask_or(circle.astype(bool), image.mask)
#         ima = np.ma.masked_array(image, mask=mask_prod)
        ima = _double(image, radius_m, x[p], y[p])
        valid_point = ima.compressed().shape[0]
        if valid_point >= th_validPoint:
            rms[i] = np.std(ima)
        else:
            pass
    ord_rms = np.sort(rms)
    return rms, ord_rms

def _double(image, radius_m, x, y):
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
    Returns
    -------
        final_ima: numpy masked array
            circle cropped image
    '''
    ps = 1/OttParameters.pscale
    if radius_m<=0.1:
        r1 = 0.1
        r_px1 = r1/ps
        ima = _circleImage(image, x, y, r_px1)
        new_ima = surf_fit(ima)
        r2 = radius_m
        r_px2 = r2/ps
        final_ima = _circleImage(new_ima, x, y, r_px2)
    else:
        r1 = radius_m
        r_px1 = r1/ps
        ima = _circleImage(image, x, y, r_px1)
        final_ima = surf_fit(ima)
    return final_ima

def _circleImage(image, x, y, raggio_px):
    circle = np.ones((image.shape[0],image.shape[1])).astype(int)
    rr, cc = draw_circle(x, y, raggio_px)
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
    zOnM4 = ZernikeOnM4()
    coef, mat = zOnM4.zernikeFit(ima,
                                np.array([2, 3]))
    new_image = zOnM4.zernikeSurface(coef, ima.mask, mat)
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
    ax = ((image - shift(image, [1,0]))/OttParameters.pscale)
    ay = ((image - shift(image, [0,1]))/OttParameters.pscale)
    mask = np.ma.mask_or(image.mask, shift(image.mask, [1,0]))
    sp = np.sqrt(((ax))**2+((ay))**2)
    #s = np.sqrt((np.arctan(ax))**2+(np.arctan(ay))**2)
    masked_slope = np.ma.masked_array(sp, mask=mask)
    return masked_slope


### REQ ###

def test242(image):
    sp = slope(image)
    # sp in pixel * fattore di conversione da rad ad arcsec
    slope_arcsec = sp * 1e-3 * 206265
    rms = slope_arcsec.std()
    return rms

def test243(image):
    diameter = 0 #average inter-actuator spacing 
    rms, rms_ord = patches_extractor(image, diameter)
    return rms

def test283(image):
    diameter = 0.08
    alpha, beta = curv_fit(image, diameter)
    raggio = roc(diameter, alpha, beta)
    return raggio

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
    image = np.zeros((dim, dim))
    xx = np.reshape(ogg[:,1], [dim,dim])
    yy = np.reshape(ogg[:,0], [dim,dim])
    z = ogg[:,2]


