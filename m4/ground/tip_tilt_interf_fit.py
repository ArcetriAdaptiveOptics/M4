
import numpy as np

def fit(image):
    dim_x, dim_y = [image.shape[0], image.shape[1]]
    mask = np.invert(image.mask).astype(int)
    mm = image.mask.astype(int)
    ll = np.array(np.where(mm[250]==0))
    ll1 = ll[0]
    pixel_dim = ll1[ll1.size-1] - ll1[0]

    x = np.linspace(1, dim_x, dim_x)
    y =  np.linspace(1, dim_y, dim_y)
    xv, yv = np.meshgrid(y, x)

    aa = xv[np.where(mask==1)]
    bb = yv[np.where(mask==1)]
    cc = np.ones(aa.size)
    z = np.zeros((aa.size, 3))
    z[:,0] = aa
    z[:,1] = bb
    z[:,2] = cc

    inv = np.linalg.pinv(z)
    c = np.dot(inv, image.compressed())
    #c slope per pixel
    #c*pixel/4/lamba*2

    coef = np.array([c[0], c[1]])
    interf_coef = coef * pixel_dim / 2 / (633e-9)
    return coef, interf_coef

# file_name = '/home/labot/immagine.fits'