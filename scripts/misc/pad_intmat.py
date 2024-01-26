'''
Authors
  - R. Briguglio: written in 2022
'''

from m4.ground import read_data
import numpy as np
from m4.ground import zernike as zern
from m4.ground import geo as geo
from astropy.io import fits
ic=read_data.InterferometerConverter()
base = '/mnt/m4storage/Data/M4Data/OPTData/PARTest/20220802/'
#base = '/mnt/data/M4/Data/M4Data/OPTData/PARTest/20220802/'


def parzmodes(fname):
    base = '/mnt/m4storage/Data/M4Data/OPTData/PARTest/20220802/'
    img=ic.fromPhaseCam6110(base+fname)
    #x0,y0,r,xx,yy = geo.qpupil(img.mask)
    mask = np.invert(img.mask).astype(int)
    #cx,cy = [950,950]
    #mask1 = geo.draw_mask(mask, cx,cy,240, out=1)
    #img2= np.ma.masked_array(img.data, -1*mask1+1)
    zlist = [1,2,3,4,5,6,7,8,9,10,11]
    #ccx, zmatx = zern.zernikeFit(img2, zlist)
    ccx, zmatx = zern.zernikeFit(img, zlist)
    #clf(); imshow(img); colorbar(); title(fname)
    ast = ccx[4:6]
    tref=ccx[8:10]
    cc = np.array((ast, tref))
    return cc

#names convention
# p0   ref for pad2
#p1,p5, pad2, X+, Y+
#p0bis  ref for pad 1,3
#p6,p7, pad1, X+, Y+
#p8,p9, pad3, X+, Y+

#pads and DoF orientation
#         ------------    PAR cover
#              PAD2         ^ +Y    > +X  orientations are preserved (pad to pad) clock-wise. +Y is always out-ward, +X is tangential, right-ward)


#        PAD3       PAD1     
#    .
#     \ +X  /+Y      / +X \+Y  (dot indicates the arrow, or positive sign)
#          '        '      '  


#************  PAD IFF MEASUREMENTS

p00 = ['20220802_144100.4D','20220802_144200.4D','20220802_144300.4D','20220802_144400.4D','20220802_144500.4D']  #ref for pad2
p10 = ['20220802_145900.4D','20220802_150000.4D','20220802_150100.4D','20220802_150200.4D','20220802_150300.4D'] #pad2, +x
p50 = ['20220802_154200.4D','20220802_154300.4D','20220802_154400.4D','20220802_154500.4D','20220802_154600.4D'] # pad2, +y

p0bis = ['20220804_102400.4D','20220804_102500.4D','20220804_102600.4D','20220804_102700.4D','20220804_102800.4D' ]  #ref for pad 1,3
p60 = ['20220804_105900.4D', '20220804_110000.4D', '20220804_110100.4D', '20220804_110200.4D', '20220804_110300.4D'] #pad1, +X
p70 = ['20220804_103700.4D', '20220804_103800.4D', '20220804_103900.4D', '20220804_104000.4D', '20220804_104100.4D']  #pad1, +Y
p80 = ['20220804_113900.4D','20220804_114000.4D', '20220804_114100.4D', '20220804_114200.4D', '20220804_114300.4D'] #pad3, +X
p90 = ['20220804_113000.4D', '20220804_113100.4D', '20220804_113200.4D', '20220804_113300.4D', '20220804_113400.4D'] #pad3, +Y



#*********** COMPUTING ZERNIKEs
pz00 = np.zeros((4,5)) #nmodes_nmeas
pz10 = np.zeros((4,5))
pz50 = np.zeros((4,5))
pz0bis = np.zeros((4,5))
pz60 = np.zeros((4,5))
pz70 = np.zeros((4,5))
pz80 = np.zeros((4,5))
pz90 = np.zeros((4,5))
z2flat = np.zeros((4,5))

for ii in range(5):
    pz00[:,ii] = parzmodes(p00[ii]).flatten()
    pz10[:,ii] = parzmodes(p10[ii]).flatten()
    pz50[:,ii] = parzmodes(p50[ii]).flatten()
    pz0bis[:,ii] = parzmodes(p0bis[ii]).flatten()    
    pz60[:,ii] = parzmodes(p60[ii]).flatten()
    pz70[:,ii] = parzmodes(p70[ii]).flatten()
    pz80[:,ii] = parzmodes(p80[ii]).flatten()
    pz90[:,ii] = parzmodes(p90[ii]).flatten()

pz0m = np.mean(pz00, 1)
pz1m = np.mean(pz10, 1)
pz5m = np.mean(pz50, 1)
pz0bism = np.mean(pz0bis,1)
pz6m = np.mean(pz60, 1)
pz7m = np.mean(pz70, 1)
pz8m = np.mean(pz80, 1)
pz9m = np.mean(pz90, 1)


z1m = pz1m - pz0m
z5m = pz5m - pz0m
z6m = pz6m -pz0bism
z7m = pz7m -pz0bism
z8m = pz8m -pz0bism
z9m = pz9m -pz0bism

#**** this is the intmat
cmat = np.zeros((6,4))
cmat[0,:] = z1m
cmat[1,:] = z5m
cmat[2,:] = z6m
cmat[3,:] = z7m
cmat[4,:] = z8m
cmat[5,:] = z9m
cmat = cmat.T


invmat = np.linalg.pinv(cmat)

