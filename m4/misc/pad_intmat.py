'''
Autore: Runa
matrice interazione pad per correggere astigmatismo parabola
'''
import numpy as np
from matplotlib import pyplot as plt
from m4.ground import read_data
from m4.ground import zernike as zern
from m4.ground import geo as geo
from astropy.io import fits


ic=read_data.InterferometerConverter()
base = '/mnt/data/M4/Data/M4Data/OPTData/PAR_Meas/20220802/'

def parzmodes(fname):
#    base = '/mnt/data/M4/Data/M4Data/OPTData/PAR_Meas/20220620/'
    img=ic.fromNew4D(base+fname)
    #x0,y0,r,xx,yy = geo.qpupil(img.mask)
    mask = np.invert(img.mask).astype(int)
    cx,cy = [950,950]
    mask1 = geo.draw_mask(mask, cx,cy,240, out=1)
    img2= np.ma.masked_array(img.data, -1*mask1+1)
    zlist = [1,2,3,4,5,6,7,8,9,10,11]
    ccx, zmatx = zern.zernikeFit(img2, zlist)
    plt.clf()
    plt.imshow(img2)
    plt.colorbar()
    plt.title(fname)
    foc = ccx[3]
    ast = ccx[4:6]
    tref=ccx[8:10]
    cc = np.array((foc, ast, tref))
    return cc

def rotmat(beta):
    theta = beta*np.pi/180
    c, s = np.cos(theta), np.sin(theta)
    rr = np.array(((c, -s), (s, c)))
    return rr

fname = '20220802_144100.4D'
# Linearity test
P0 = '20220802_144100.4D' # 0,0
P1 = '20220802_145900.4D' # X
P5 = '20220802_154200.4D' # Y

p00 = ['20220802_144100.4D','20220802_144200.4D','20220802_144300.4D','20220802_144400.4D','20220802_144500.4D']
p10 = ['20220802_145900.4D','20220802_150000.4D','20220802_150100.4D','20220802_150200.4D','20220802_150300.4D']
p50 = ['20220802_154200.4D','20220802_154300.4D','20220802_154400.4D','20220802_154500.4D','20220802_154600.4D']

p0bis = []
p60 = []
p70 = []


pz00 = np.zeros((5,5)) #nmodes_nmeas
pz10 = np.zeros((5,5))
pz50 = np.zeros((5,5))
pz0bis = np.zeros((5,5))
pz60 = np.zeros((5,5))
pz70 = np.zeros((5,5))

for ii in range(5):
    pz00[:,ii] = parzmodes(p00[ii]).flatten()
    pz10[:,ii] = parzmodes(p10[ii]).flatten()
    pz50[:,ii] = parzmodes(p50[ii]).flatten()
    pz0bis[:,ii] = parzmodes(p0bis[ii]).flatten()
    pz60[:,ii] = parzmodes(p60[ii]).flatten()
    pz70[:,ii] = parzmodes(p70[ii]).flatten()

pz0m = np.mean(pz00, 1)
pz1m = np.mean(pz10, 1)
pz5m = np.mean(pz50, 1)
pz0bism = np.mean(pz0bis,1)
pz6m = np.mean(pz60, 1)
pz7m = np.mean(pz70, 1)

z1m = pz1m - pz0m
z5m = pz5m - pz0m
z6m = pz6m -pz0bism
z7m = pz7m -pz0bism

cmat = np.zeros((4,5))
cmat[0,:] = z1m
cmat[1,:] = z5m
cmat[2,:] = z6m
cmat[3,:] = z7m
cmat = cmat[:,1::] #delete this line for focus
cmat = cmat.T

invmat = np.linalg.pinv(cmat)
cmd0 = np.matmul(invmat, pz0m)
pz0m - np.matmul(cmat, cmd0)
cmd1 = np.array((cmd0, cmd0))

rmat = rotmat(120)
cmat1 = np.matmul(cmat, rmat.T)

# ccm = np.concatenate((cmat, mat1), 1)
#  
# invmat2 = np.linalg.pinv(ccm)
# cmd2 = np.matmul(invmat2, pz0m)






