#script for the comparison of the PAR shape before-after coating
from m4.ground import zernike as zern
from m4.ground import geo as geo
from matplotlib.pyplot import *
from astropy.io import fits
from m4.ground import timestamp
from m4.ground import read_data
from m4.mini_OTT import timehistory as th
import numpy as np
import os



tns = '20230127_163227'
tns = '20230401_134742'
nf = 200
fl = th.fileList(tns)
img0 = th.averageFrames(0,nf,fl)
img0 = th.removeZernike(img0,[1,2,3,4,7,8])
print(img0.std())

img1 = th.averageFrames(nf+1,nf+nf,fl)
img1 = th.removeZernike(img1,[1,2,3,4,7,8])
dd = img1-img0
print(dd.std())

