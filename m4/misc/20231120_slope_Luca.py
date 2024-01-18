from m4.configuration import start
from m4.mini_OTT import timehistory as th
from m4.ground import zernike
from m4.ground import timestamp
import time as tt
from m4 import noise
import os
import numpy as np
from matplotlib.pyplot import *
from m4.ground import geo
from importlib import reload
from m4.misc import image_registration_lib as imgreg
from astropy.io import fits as pyfits
import scipy


def congrid2D(img, dimx,dimy, method='linear', centre=False, minusone=False):
    dims = img.shape
    xo = np.linspace(1,dims[0], dims[0])
    yo = np.linspace(1,dims[1], dims[1])
    x = np.linspace(1,dims[0],int(dimx))
    y = np.linspace(1,dims[1],int(dimy))
    xg,yg =np.meshgrid(x,y,indexing='ij')
    interp = scipy.interpolate.RegularGridInterpolator((xo,yo),img,method='linear')
    newimg = interp((xg,yg))
    return newimg

def mask_linear(im,x1, y1, x2, y2,d):

    # d=10
    x = np.arange(1,np.shape(im)[0]+1,1)
    y = np.arange(1,np.shape(im)[1]+1,1)
    xx, yy = np.meshgrid(x, y)
    m=(y2-y1)/(x2-x1)
    #m=(x2-x1)/(y2-y1)
    q=y1-m*x1
    
    temp1=yy-m*xx-q-d/2>0
    temp2=yy-m*xx-q+d/2<0
    mask=np.array(temp1+temp2)
    #plt.figure(); plt.imshow(mask); plt.show()
    
    return mask

def circ_mask(im,x0,y0,R,maschero=1):
    
    mm=im.mask
    
    # definisco la maschera ma
    nx=np.shape(im)[0]
    ny=np.shape(im)[1]
    x = np.linspace(1, nx, nx)
    y = np.linspace(1, ny, ny)
    xv, yv = np.meshgrid(x, y)
    if maschero==1:
        temp=(xv-x0)**2+(yv-y0)**2<R**2
    else:
        temp=(xv-x0)**2+(yv-y0)**2>R**2
                
    
    # figure(); imshow(temp); colorbar(); show()
    mask=np.logical_or(mm,temp)
    # figure(); imshow(mask); colorbar(); show()
    return mask

# def compSlope(img,px, rfact):
#     ss=np.array(np.shape(img))
#     sli = img-np.roll(img,(1,1),axis=(0,1))/px
#     slid = geo.congrid2D(sli,(ss/rfact).astype(int))
#     slim = geo.congrid2D(-1*sli.mask+1,(ss/rfact).astype(int))
#     sli = np.ma.masked_array(slid,(-1*slim+1))
#     return sli.std()

##
dove="/home/m4/Downloads"
name="20231029_150436-OTT-Cal.fits"
fits_file_name = os.path.join(dove, name)
hduList=pyfits.open(fits_file_name)
image = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))

## rebin -> slope, con rebin di Xompero
close('all')
im0=image.copy()
im=im0

figure(); imshow(im); colorbar(); show()


k=4
px=0.0007
rfact=k

bin=rebin(im, (np.shape(im)[0]/k,np.shape(im)[1]/k) )
bin_mask=rebin(im.mask, (np.shape(im)[0]/k,np.shape(im)[1]/k)  )
sli = (bin-np.roll(bin,(1,1),axis=(0,1)))/(px*rfact)
sli.mask=bin_mask
# 

mm=circ_mask(sli,1210/k,766/k,4); sli.mask=mm 
mm=circ_mask(sli,820/k,760/k,4); sli.mask=mm
mm=circ_mask(sli,1020/k,434/k,4); sli.mask=mm
mm=circ_mask(sli,1020/k,506/k,350/k,0); sli.mask=mm
rms1=np.std(sli)
media=np.average(sli)

print('rebin -> slope:', rms1*180/np.pi*3600)
figure(); imshow(sli); colorbar(); show()

# slope -> rebin
im=image.copy()
sli = (im-np.roll(im,(1,1),axis=(0,1)))/(px)
bin=rebin(sli,(np.shape(im)[0]/k,np.shape(im)[1]/k) )
bin.mask=rebin(im.mask, (np.shape(im)[0]/k,np.shape(im)[1]/k) )


mm=circ_mask(bin,1210/k,766/k,4); bin.mask=mm 
mm=circ_mask(bin,820/k,760/k,4); bin.mask=mm
mm=circ_mask(bin,1020/k,434/k,4); bin.mask=mm
mm=circ_mask(bin,1020/k,506/k,350/k,0); bin.mask=mm

rms2=np.std(bin)

figure(); imshow(bin); show()

print('slope -> rebin:', rms2*180/np.pi*3600)


## rebin -> slope, con rebin di runa
# im=image.copy()
# k=4
# px=0.0007
# rfact=k
# 
# bin=congrid2D(im,2048/k,2048/k)
# sli = (bin-np.roll(bin,(1,1),axis=(0,1)))/(px*rfact)
# rms=np.std(sli)
# 
# # slope -> rebin
# im=image.copy()
# sli = (im-np.roll(im,(1,1),axis=(0,1)))/px
# bin=congrid2D(sli,2048/k,2048/k)
# rms2=np.std(bin)
# 
# 
# print('rebin -> slope:',rms*180/np.pi*3600)
# print('slope -> rebin:',rms2*180/np.pi*3600)



##
def rebin(a, new_shape, sample=False):
    """
    Replacement of IDL's rebin() function for 2d arrays.

    Resizes a 2d array by averaging or repeating elements.
    New dimensions must be integral factors of original dimensions,
    otherwise a ValueError exception will be raised.

    Parameters
    ----------
    a : ndarray
        Input array.
    new_shape : 2-elements sequence
        Shape of the output array
    sample : bool
        if True, when reducing the array side elements are set
        using a nearest-neighbor algorithm instead of averaging.
        This parameter has no effect when enlarging the array.

    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array  the data are averaged,
        unless the sample parameter is set.
        If the new shape is bigger array elements are repeated.

    Raises
    ------
    ValueError
        in the following cases:
         - new_shape is not a sequence of 2 values that can be converted to int
         - new dimensions are not an integral factor of original dimensions
    NotImplementedError
         - one dimension requires an upsampling while the other requires
           a downsampling

    Examples
    --------
    >>> a = np.array([[0, 1], [2, 3]])
    >>> b = rebin(a, (4, 6)) #upsize
    >>> b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])
    >>> rebin(b, (2, 3)) #downsize
    array([[0. , 0.5, 1. ],
           [2. , 2.5, 3. ]])

    >>> rebin(b, (2, 3)) #downsize
    array([[0. , 0.5, 1. ],
           [2. , 2.5, 3. ]])
    >>> rebin(b, (2, 3), sample=True) #downsize
    array([[0, 0, 1],
           [2, 2, 3]])
    """

    # unpack early to allow any 2-length type for new_shape
    m, n = map(int, new_shape)

    if a.shape == (m, n):
        return a

    M, N = a.shape

    if m <= M and n <= M:
        if (M//m != M/m) or (N//n != N/n):
            raise ValueError('Cannot downsample by non-integer factors')

    elif M <= m and M <= m:
        if (m//M != m/M) or (n//N != n/N):
            raise ValueError('Cannot upsample by non-integer factors')

    else:
        raise NotImplementedError('Up- and down-sampling in different axes '
                                  'is not supported')

    if sample:
        slices = [slice(0, old, float(old) / new)
                  for old, new in zip(a.shape, (m, n))]
        idx = np.mgrid[slices].astype(np.int32)
        return a[tuple(idx)]
    else:
        if m <= M and n <= N:
            return a.reshape((m, M//m, n, N//n)).mean(3).mean(1)
        elif m >= M and n >= M:
            return np.repeat(np.repeat(a, m/M, axis=0), n/N, axis=1)