'''
This code is meant to perform other analysis on the noise data processed with noise.noise_vibrations and noise.convection_noise.

author
L. Oggioni

'''

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
from m4.utils import image_registration_lib as imgreg
from astropy.io import fits as pyfits

def _funFit(x, a, b, c):
    fun = -np.exp(-a * (x - b)) + c
    return fun
    
def _curvFit(param, x, rms_nm):
    from scipy.optimize import curve_fit
    pp, pcov = curve_fit(_funFit, x, rms_nm, param)
    fit = _funFit(x, *pp)
    return pp, fit

##
path='/mnt/m4storage/Data/Results/Noise'

tn=[]
tn.append('20231011_144711')
tn.append('20231007_113703')
tn.append('20231006_232601')
tn.append('20231028_084255')
##
close('all')
tt=[]; RMS=[]; xx=[]
for j in np.arange(0,np.size(tn),1):
#
    file=os.path.join(path,tn[j],'tiptilt_vector_conv.fits')
    hduList = pyfits.open(file)
    tt.append(hduList[0].data)
    
    file=os.path.join(path,tn[j],'rms_vector_conv.fits')
    hduList = pyfits.open(file)
    RMS.append(hduList[0].data)
    
    file=os.path.join(path,tn[j],'time_vector_conv.fits')
    hduList = pyfits.open(file)
    xx.append(hduList[0].data)


    figure(1); plot(xx[j], RMS[j],'-o',label=tn[j]); legend()
    xlabel('time [s]'); ylabel('rms [nm]'); title('structure function')
    
    x=xx[j].copy()
    rms=RMS[j].copy()
    rms_nm = rms * 1e9
    param = [5, 0.5, 32]
    try:

        pp, fit = _curvFit(param, x, rms_nm)
        decorr_time = 1 / pp[0] + pp[1]
    except:
        pp = np.array([0, 0, rms[-1] * 1e9])
        decorr_time = -1
        fit = rms_nm.copy() * 0
    figure()
    clf()
    plot(x, rms * 1e9, '-o', label='meas')
    xlabel('time [s]')
    ylabel('rms [nm]')
    plot(x, fit, '-', label='fit')
    grid()
    plot([x[0], x[-1]], [pp[2], pp[2]], '--r', linewidth=3,label='%.2f [nm]' % pp[2])
    legend()
    title(tn[j]+' structure function')

show()
##
close('all')
tt=[]; RMS=[]; n_temp=[]
for j in np.arange(0,np.size(tn),1):

    file=os.path.join(path,tn[j],'tiptilt_vector_0.fits')
    hduList = pyfits.open(file)
    tt.append(hduList[0].data)
    
    file=os.path.join(path,tn[j],'rms_vector_0.fits')
    hduList = pyfits.open(file)
    RMS.append(hduList[0].data)
    
    file=os.path.join(path,tn[j],'n_temp_vector_0.fits')
    hduList = pyfits.open(file)
    n_temp.append(hduList[0].data)

    figure(1); plot(n_temp[j], RMS[j] * 1e9,'-o',label=tn[j]); legend()
    xlabel('n_temp'); ylabel('rms [nm]')

    figure(2); plot(n_temp[j], tt[j] * 1e9,'-o',label=tn[j]); legend()
    xlabel('n_temp'); ylabel('tt [nm]')

    figure()
    plot(n_temp[j], RMS[j] * 1e9, '-o',label='rms [nm]')
    plot(n_temp[j], tt[j] * 1e9, '-o',label='TipTilt [nm]')
    xlabel('n_temp')
    ylabel('[nm]')
    title(tn[j])
    legend()
    grid()

show()
    
