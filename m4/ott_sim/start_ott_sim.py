# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:06:35 2020

@author: Runa
list of startup instruction to startup the simulated OTT and get interferometer frames
usage:
    import everything as follows
    setup the OTT:
        ott.slide, ott.rslide, ott.angle
            slide=0.8 for segment, rslide =0.3 for alignment, angle 0,60,120... for segments
    setup the mirrors in ott:
        ott.param, ott.m4, ott.refflat

    acquire an image
"""

import numpy as np
from m4.ground import geo

### CHIARA ###
from m4.utils.interface_4D import comm4d
c4d = comm4d()
from m4.configuration import start
ott=start.create_ott()

ott.slide(0.75)
ott.angle(70.)
ott.rslide(0.3)
ott.parab([0,0,1e-5,0,0,0])
ott.m4([0,0,0,1e-6,0,0])
ott.refflat([0,0,0,0,1e-6,0])

p,m = c4d.acq4d(ott, 1, show=1)

from m4.ott_sim.ott_images import OttImages
ottIma = OttImages(ott)

p,m=ottIma.ott_smap(show=1)


a=np.zeros((21,21)) 
b=geo.draw_mask(a,10,10,10)

segmask1 = np.ma.make_mask(ott.segmask1)



'''
### RUNA ###
from importlib import reload
import os
from matplotlib import pyplot as plt
a='D:\Astro\ARCETRI\Python\M4-master'
os.chdir(a)
from m4.configuration.create_ott import *
from m4.configuration import start
ott=start.create_ott()
from m4.ott_sim.ott_images import *
from m4.utils.interface_4D import comm4d
#test
plt.imshow(ott.mask)



ott.slide(0.75)
ott.angle(30.)
ott.rslide(0.3)
ott.parab([0,0,1e-5,0,0,0])
ott.m4([0,0,0,1e-6,0,0])
ott.refflat([0,0,0,0,1e-6,0])


p,m = comm4d.acq4d(1,show=1)
#which is a call to the simulated interf.acquisition
p,m=ott_smap(show=1)



a=np.array([20,20])
b=geo.draw_mask(a,10,10,10)
plt.imshow(b)


bname =conf.path_name.MIRROR_FOLDER+'/'+conf.mirror_conf
fname = bname+'/py-sect4-mask.fits'
hduList = pyfits.open(fname)
m = hduList[0].data
segmask1 = ma.make_mask(m)

bname =conf.path_name.MIRROR_FOLDER+'/'+conf.mirror_conf
fname = bname+'/if_sect4_rot-bin2.fits'
hduList = pyfits.open(fname)
ifmat = hduList[0].data

 s1 = segmask1.copy()
 s1=s1.astype(float)
 b[segmask1 == True]=ifmat[3,:]  #if zonale

 comm  = np.dot(ott.ifmat.T,ott.vmat[:,3]) #ifmodale
 s1 = segmask1.copy()
 s1=s1.astype(float)
 s1[segmask1 == True]=comm
 plt.imshow(s1)
'''
