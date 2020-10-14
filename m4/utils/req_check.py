'''
@author: rm
'''

from m4.configuration.ott_parameters import *
from astropy.io import fits as pyfits
from skimage.draw import circle as draw_circle
from m4.utils.zernike_on_m_4 import ZernikeOnM4
from scipy.ndimage.interpolation import shift

def patches_extractor(image, metri):
    #metri = 0.05
    ps = 1/OttParameters.PIXEL_SCALE
    raggio = metri/ps

    nn = image.compressed().shape[0]
    th_validPoint = np.pi * raggio**2 * 30 / 100
    idx = np.where(image.mask==0)
    x = idx[0]
    y = idx[1]

    cube = None
    for i in range(200):
        print('%d' %i)
        circle = np.ones((image.shape[0],image.shape[1])).astype(int)
        rr, cc = draw_circle(x[i], y[i], raggio)
        circle[rr, cc] = 0
        mask_prod = np.ma.mask_or(circle.astype(bool), image.mask)
        ima = np.ma.masked_array(image, mask=mask_prod)
        valid_point = ima.compressed().shape[0]
        if valid_point >= th_validPoint:
            if cube is None:
                cube = ima
            else:
                cube = np.ma.dstack((cube, ima))
        else:
            pass
    return cube

def surf_fit(cube):
    zOnM4 = ZernikeOnM4()
    cube_tt = None
    for i in range(cube.shape[2]):
        print('%d' %i)
        image = cube[:,:,i]
        coef, mat = zOnM4.zernikeFit(image,
                                    np.array([2, 3]))
        new_image = zOnM4.zernikeSurface(coef, image.mask, mat)
        if cube_tt is None:
            cube_tt = new_image
        else:
            cube_tt = np.ma.dstack((cube_tt, new_image))
    return cube_tt

def std_calculator(cube):
    tt = np.zeros(cube.shape[2])
    for i in range(cube.shape[2]):
        tt[i] = np.std(cube[:,:,i])
    ord_tt = np.sort(tt)
    return tt, ord_tt

def curv_fit(image, test_diameter):
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
    wfe = test_diameter**2 /(8*np.sqrt(3)) * np.sqrt(2*(alpha-beta)**2 - (alpha+beta)**2)
    rho = test_diameter/2
    raggio = (rho**2 + wfe**2)/(2*wfe)
    return raggio

def slope(image):
    ax = ((image - shift(image, [1,0]))/OttParameters.PIXEL_SCALE)
    ay = ((image - shift(image, [0,1]))/OttParameters.PIXEL_SCALE)
    sp = np.sqrt(((ax))**2+((ay))**2)
    s = np.sqrt((np.arctan(ax))**2+(np.arctan(ay))**2)
    return s

### REQ ###

def test242(image):
    s = slope(image)
    #C = (dimensione pixel[mm] x 206265) / Focale utilizzata
    slope_arcsec = s * 1e-3 * 206265
    rms = slope_arcsec.std()
    return rms

def test243(image):
    diameter = 0 #average inter-actuator spacing 
    cube = patches_extractor(image, diameter)
    rms_vector = np.zeros(cube.shape[2])
    for i in range (cube.shape[2]):
        ima = cube[:,:,i]
        rms_vector[i] = ima.std()
    return rms_vector

def test283(image):
    test_diameter = 0.08
    alpha, beta = curv_fit(image, test_diameter)
    raggio = roc(test_diameter, alpha, beta)
    return raggio

def imaTest():
    ff = '/Users/rm/Desktop/Arcetri/M4/ProvaCodice/ZernikeCommandTest/20191210_110019/mode0002_measure_segment00_neg.fits'
    hduList=pyfits.open(ff)
    image = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
    return image