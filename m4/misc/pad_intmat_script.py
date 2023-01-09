from m4.ground import read_data
from m4.ground import zernike as zern
from m4.ground import geo as geo
from astropy.io import fits
from m4.misc import pad_intmat as padim

ic=read_data.InterferometerConverter()
base = '/mnt/m4storage/Data/M4Data/OPTData/PARTest/'
#base = '/mnt/data/M4/Data/M4Data/OPTData/PARTest/20220802/'

#restoring the pad intmat
im = padim.invmat

def getparzern(tn):
    fold = tn[0:8]
    img = read_data.read_phasemap(base+fold+'/'+tn+'.fits')
    zlist = [1,2,3,4,5,6,7,8,9,10,11]
    #ccx, zmatx = zern.zernikeFit(img2, zlist)
    ccx, zmatx = zern.zernikeFit(img, zlist)
    print(ccx)
    ast = ccx[4:6]
    return ast





def rms(val):
    t=mean(val**2)
    return sqrt(t)

def savenp(fname):
    basein = '/mnt/m4storage/Data/M4Data/OPTData/PARTest/20220802/'
    baseout ='/home/m4/Desktop/data4ADS/PARTest/20220802/'
    img=ic.fromPhaseCam6110(base+fname)
    np.save(baseout+fname,img.data)

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
    clf(); imshow(img); colorbar(); title(fname)
    ast = ccx[4:6]
    tref=ccx[8:10]
    cc = array((ast, tref))
    return cc

def rotmat(beta):
    theta = beta*np.pi/180
    c, s = np.cos(theta), np.sin(theta)
    rr = np.array(((c, -s), (s, c)))
    return rr

fname = '20220802_144100.4D'
# Linearity test
P0 = '20220802_144100.4D' # 0,0 , ref for pad2
P1 = '20220802_145900.4D' # X, pad2
P5 = '20220802_154200.4D' # Y, pad2

#intmat test
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


#************ image to be FLATTENED
img2flat = ['20220804_102400.4D','20220804_102500.4D','20220804_102600.4D','20220804_102700.4D','20220804_102800.4D' ]  # corresponda to p00bis

img2flat2 = ['20220804_124500.4D', '20220804_124600.4D', '20220804_124700.4D', '20220804_124800.4D', '20220804_124900.4D']

img2flat3 = ['20220804_142800.4D', '20220804_142900.4D','20220804_143000.4D','20220804_143100.4D', '20220804_143200.4D']

img2flat4 = ['20220804_153800.4D','20220804_153900.4D', '20220804_154000.4D', '20220804_154100.4D','20220804_154200.4D']

img2flat4bis = ['20220804_155900.4D','20220804_160000.4D','20220804_160100.4D','20220804_160200.4D','20220804_160300.4D']  # no cmd applied, for repatab./stability vs imgflat4

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

pz0m = mean(pz00, 1)
pz1m = mean(pz10, 1)
pz5m = mean(pz50, 1)
pz0bism = mean(pz0bis,1)
pz6m = mean(pz60, 1)
pz7m = mean(pz70, 1)
pz8m = mean(pz80, 1)
pz9m = mean(pz90, 1)

#*** Testing the meas REPEATABILITY
pz0s = std(pz00, 1)
pz1s = std(pz10, 1)
pz5s = std(pz50, 1)
pz0biss = std(pz0bis,1)
pz6s = std(pz60, 1)
pz7s = std(pz70, 1)
pz8s = std(pz80, 1)
pz9s = std(pz90, 1)


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


#****** computing the ZERNIKE to flatten
z2flat = np.zeros((4,nimg))
nimg = len(img2flat)
for ii in range(nimg):
    z2flat[:,ii] = parzmodes(img2flat[ii]).flatten()
pz2flat = mean(z2flat, 1)

nimg = len(img2flat2)
z2flat2 = np.zeros((4,nimg))
for ii in range(nimg):
    z2flat2[:,ii] = parzmodes(img2flat2[ii]).flatten()
pz2flat2 = mean(z2flat2, 1)

nimg = len(img2flat3)
z2flat3 = np.zeros((4,nimg))
for ii in range(nimg):
    z2flat3[:,ii] = parzmodes(img2flat3[ii]).flatten()
pz2flat3 = mean(z2flat3, 1)

nimg = len(img2flat4)
z2flat4 = np.zeros((4,nimg))
for ii in range(nimg):
    z2flat4[:,ii] = parzmodes(img2flat4[ii]).flatten()
pz2flat4 = mean(z2flat4, 1)

nimg = len(img2flat4bis)
z2flat4bis = np.zeros((4,nimg))
for ii in range(nimg):
    z2flat4bis[:,ii] = parzmodes(img2flat4bis[ii]).flatten()
pz2flat4bis = mean(z2flat4bis, 1)

nimg = len(img2flat5)
z2flat5 = np.zeros((4,nimg))
for ii in range(nimg):
    z2flat5[:,ii] = parzmodes(img2flat5[ii]).flatten()
pz2flat5 = mean(z2flat5, 1)



#here the flattening command is computed
invmat = np.linalg.pinv(cmat)
cmd0 = -matmul(invmat, pz2flat)
cmd1 = -matmul(invmat, pz2flat2)
cmd2 = -matmul(invmat, pz2flat3)
cmd3 = -matmul(invmat, pz2flat4)
cmd5 = -matmul(invmat, pz2flat5)

#******* evaluating the expected residual
pz2flat -matmul(cmat, cmd0)


#********* saving the data for Matteo
img2flat = ['20220804_102400.4D','20220804_102500.4D','20220804_102600.4D','20220804_102700.4D','20220804_102800.4D' ]  # corresponda to p00bis

img2flat2 = ['20220804_124500.4D', '20220804_124600.4D', '20220804_124700.4D', '20220804_124800.4D', '20220804_124900.4D']

img2flat3 = ['20220804_142800.4D', '20220804_142900.4D','20220804_143000.4D','20220804_143100.4D', '20220804_143200.4D']

img2flat4 = ['20220804_153800.4D','20220804_153900.4D', '20220804_154000.4D', '20220804_154100.4D','20220804_154200.4D']

img2flat4bis = ['20220804_155900.4D','20220804_160000.4D','20220804_160100.4D','20220804_160200.4D','20220804_160300.4D']

img2flat5 = ['20221216_104500-16ave.4D','20221216_104000-16ave.4D' ]  # corresponda to p00bis

#for ii in range(5):
#    savenp(img2flat[ii])
#    savenp(img2flat2[ii])
#    savenp(img2flat3[ii])
#    savenp(img2flat4[ii])
#    savenp(img2flat4bis[ii])


