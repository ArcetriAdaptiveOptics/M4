from m4.configuration import start
from m4.ground import zernike as zern
from m4.ground import geo as geo
from matplotlib.pyplot import *
from astropy.io import fits
from m4.ground import timestamp 
from m4.ground import read_data
import numpy as np
import os
#ic=read_data.InterferometerConverter()

conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
ott, interf, dm = start.create_ott(conf)
ts = timestamp.Timestamp()

#salva sul computer locale
base = '/mnt/m4storage/Data/M4Data/OPTData/PARTest/'

def acq_par(nframes=16, delay=1):
    thefold = ts.today()
    myfold=base+thefold
    if not os.path.exists(myfold):
        os.mkdir(myfold)
    tn = ts.now()+'.fits'
    img = interf.acquire_phasemap(nframes, delay)
    interf.save_phasemap(myfold,tn,img)  
    print(myfold)
    print(tn)    

    return img,  tn

def get_zern(im):
    #myfold =base+tn[0:8]
    #im=read_data.read_phasemap(os.path.join(myfold,tn))
    zlist = [1,2,3,4,5,6,7,8,9,10,11]
    cc, zmat = zern.zernikeFit(im, zlist)

    #print(cc[4:6])   

    return cc

def removeZernike(ima, modes=np.array([1,2,3,4])):
    coeff, mat = zernike.zernikeFit(ima, modes)
    surf = zernike.zernikeSurface(ima, coeff, mat)
    new_ima = ima-surf
    
    return new_ima, coeff

def get_average(fold):
    where = '/mnt/m4storage/Data/M4Data/OPTData/OPD_series/'
    fname_zern = 'zernike.fits'
    file_name = where + fold + '/' + fname_zern
    hdul = fits.open(file_name)
    data = hdul[0].data
    time = data[:,0]
    coeff = data[:,1::]
    ave = np.average(coeff,0)
    sstd = np.std(coeff,0)
    figure(figsize=(8,5))
    plot(coeff[:,4],'-o', label='ast X')
    plot(coeff[:,5],'-o', label='ast Y')
    legend(); grid(); ylabel('nm RMS'); xlabel('frame #'); title(fold + ' Astigmatism plot')
    figure(figsize=(8,5))
    plot(coeff[:,6],'-o', label='coma X')
    plot(coeff[:,7],'-o', label='coma Y')
    legend(); grid(); ylabel('nm RMS'); xlabel('frame #'); title(fold + ' Coma plot')
    print('Average Astigmatism =  %.3e ±  %.3e, %.3e ± %.3e  m' %(ave[4], sstd[4], ave[5], sstd[5]))
    print('Average Coma =  %.3e ± %.3e, %.3e ± %.3e m' %(ave[6],sstd[6], ave[7], sstd[7]))
    return coeff, time
