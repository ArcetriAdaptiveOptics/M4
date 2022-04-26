##script for requirements verification data analysis

##OTT setup accuracy (RA-M4-0279, RA-M4-0288)

from M4.mOTT_analysis import timehistory as th
import pandas as pd
import numpy as np
from astropy.io import fits
from arte.utils.rebin import rebin
from scipy import interpolate
import matplotlib.pyplot as plt


class reqVerification():
    
    def __init__(self, id=''):
    
        if id=='':
            print("possible reqs:")
            print("RA-M4-0324: SPL accuracy")
            print("RA-M4-0270: Inter-actuator stroke")
            print("RA-M4-0274: Slope accuracy")
            print("RA-M4-0277: Curvature accuracy")
            print("RA-M4-0279: OT Setup accuracy")
            print("RA-M4-0282: Co-phasing accuracy")
            print("RA-M4-0267: Quasi-static Zernike offset accuracy")
            print("RA-M4-0292: WFE accuracy")
            return None
        elif id=='RA-M4-0279':
            self._time=0
            self._ConvectionEffect=0
            self._ThermalEffect=0
            self._VibrationEffect=0
            self._RelaySystem =0
            self._Parabola  =0
            self._ParabolaAccuracy=0
            self._FlatMirror=0
            self._ProcessingError  =0
            self._TNlist=['' ]
            self._OutTemp=self._fillTempFromDB()
        else:
            print("Not implemented yet")
            return None
        
        
    def _fillTempFromDB(self, fname='/home/labot/git/M4/m4/mOTT_analysis/RW_20210513113809_438286_5883_1.csv'):
        return pd.read_csv(fname, parse_dates=['Data-Ora'])
            
    def computeVals(self):
        pass
        
class removeDiffPist():
    
    def __init__(self, filename='./test.fits'):
        hdul = fits.open(filename)
        img = hdul[0].data
        plt.imshow(img)
        
        sx = img.shape[0]
        sy = img.shape[1]
        cx = 406
        cy = 406
        rmin = 100
        rmax = 200
        ncircles = 11
        nsamples = 1000
        
        #generate circles
        angles =np.array([ np.arange(nsamples)/(nsamples-1)*2*np.pi])
        rays = np.array([np.linspace(rmin, rmax, ncircles)])
        th = rebin(angles, [ncircles, nsamples])
        rr = rebin(np.transpose(rays),[ncircles, nsamples] )
        xi = rr*np.cos(th) + cx
        yi = rr*np.sin(th) + cy
        x = np.arange(sx)
        y = np.arange(sy)
        xx, yy = np.meshgrid(x,y)
        
        #roi detection
        from skimage import measure
        mask = img !=0
        roi = measure.label(mask)
        m1 = roi ==1;
        m2 = roi ==2;
        m3 = roi ==3;
        
        #fit each submask with Thin Plate Spline 
        
        #using the fitted surface ????
        
        
        zi = interpolate.griddata(np.array([xx.ravel(), yy.ravel()]).T, img.ravel(),(xi,yi), method='linear')
        for ii in range(ncircles):
            plt.plot(np.arange(nsamples), zi[ii,:])

    
def scaleZernike(coeff,r0,r1):
    pass

    