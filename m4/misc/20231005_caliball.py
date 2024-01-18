
from matplotlib.pyplot import *
from m4.ground import zernike as zern
import numpy as np
from m4.mini_OTT import timehistory as th
from m4.ground import geo
from m4.ground import zernike
from astropy.io import fits as pyfits



def reqcheck(image, radius, npatch, display=True):
    zlist = [1,2,3,4]
    from m4.analyzers import requirement_analyzer as ra
    idx = (np.where(image.mask==0))[0]
    nidx = np.size(idx)
    step = int(nidx / npatch)
    ww = th.removeZernike(image,zlist)
    cc = geo.qpupil(-1*ww.mask+1)
    ps = cc[2]/radius  #pix/m
    slope = ra.test242(ww, ps)
    diff_piston = ra.diffPiston(ww)
    #roc = ra.test283(ww, ps, step, display=False)
    rms31 = ra.test243(ww, 0.015, ps, step, n_patches=None)
    rr, list_ima, result_vect = ra.patches_analysis(ww, 0.015, ps, step, n_patches=None)

    #rms500 = ra.test243(ww, 0.25, ps, step, n_patches=None)
    if display is True:
        #print('Zernike to subtract = %i' %size(zlist))
        print('#######################################')
        print('############### RESULTS  on############')
        print('#######################################')
        print('slope = %.3f asec' %(slope))
        print('diff_piston = %.3e nm' %(diff_piston))
       # print('RoC = %.3e km' %(roc/1000))
        print('WFE 31 = %.3e nm' %(rms31))
        #print('WFE 500 = %.3e nm' %(rms500))

    #return slope, diff_piston, roc, rms31, rms500
    return slope, diff_piston, rms31, list_ima


tn0 = '20231004_113904' # oldL3
tn = '20231004_172617'
end = 500

fl = th.fileList(tn)
q0=th.averageFrames(0,end, fl)
ss = q0.shape


###t
#
'''
# Create and save a cube 
cube = []
for ii in range(0,end):
    print(ii)
    tmp = th.read_phasemap(fl[ii])
    tmp = th.removeZernike(tmp,[1,2,3])
    cube.append(tmp)
cube = np.ma.masked_array(cube)
where = fl[0].rsplit('/',1)[0] + '/'
fname = where + 'cube.fits'
pyfits.writeto(fname, cube.data)
'''

#Load cube
where = fl[0].rsplit('/',1)[0] + '/'
fname = where + 'cube.fits'
hdul = pyfits.open(fname)
cube = hdul[0].data



######
mm = geo.draw_mask(np.ones(ss),ss[0]/2,ss[1]/2,ss[0]/2,out=1)
q0 = np.ma.masked_array(q0, mm)
q1 = th.removeZernike(q0,[1,2,3,4,7,8])
print(q1.std())

#imshow(q1); title('Old RS, RMS=52nm SfE')
#print(q1.std())

imgrad = 1.5/2
ottrad = 0.3
pixs = 0.0007 #mm/pix
mm1 = geo.draw_mask(np.ones(ss),ss[0]/2,ss[1]/2,ottrad/pixs,out=1)
q2 = np.ma.masked_array(q0, mm1)
q2=th.removeZernike(q2,[1,2,3,4,7,8])
figure()
imshow(q2, cmap='hot')
colorbar()
title('Relay system WFE RMS = %.3e m \n  Foot print diameter = %.2f m' %(q2.std(), 2*ottrad), fontsize=18)
show()
cc, mat = zern.zernikeFit(q2,np.arange(10)+1)

ottrad = 0.7
pixs = 0.000758 #mm/pix
mm1 = geo.draw_mask(np.ones(ss),ss[0]/2,ss[1]/2,ottrad/pixs,out=1)
q2 = np.ma.masked_array(q0, mm1)
q2=th.removeZernike(q2,[1,2,3,4,7,8])
figure()
imshow(q2,cmap='hot')
colorbar()
title('Relay system RMS = %.3e m \n Average of %i images' %(q2.std(), end))
show()
cc, mat = zern.zernikeFit(q2,np.arange(10)+1)


ottrad = 0.7
pixs = 0.000758 #mm/pix
mm1 = geo.draw_mask(np.ones(ss),ss[0]/2,ss[1]/2,ottrad/pixs,out=1)
q2 = np.ma.masked_array(q0, mm1)
q2=th.removeZernike(q2,[1,2,3,4,5,6,7,8])
figure()
imshow(q2,cmap='hot')
colorbar()
title('Relay system RMS = %.3e m \n Average of %i images' %(q2.std(), end))
show()
cc, mat = zern.zernikeFit(q2,np.arange(10)+1)





#slope, diff_piston, rms31, li=reqcheck(q2,0.3,10000, display=False)



assert(False)
stdvec= []
coeff = []
for ii in range(0, len(cube)):
    print(ii)
    cc = cube[ii].data 
    idm = np.isnan(cc)
    mm = np.zeros(np.shape(cc))
    mm[idm] = 1
    aa = (-1*mm+1)*(-1*mm1+1)
    aa = -1*aa+1
    img = np.ma.masked_array(cc, aa)
    cf, mat = zern.zernikeFit(img,[1,2,3])
    coeff.append(cf)
    tmp = np.ma.masked_array(cc, mm)
    stdvec.append(tmp.std())

coeff = np.array(coeff)
sigma = np.mean(stdvec) + np.std(stdvec)
figure()
plot(np.sqrt(coeff[:,1]**2 + coeff[:,2]**2),'-o')

#plot(stdvec,'-o')
#plot([0,end],[sigma, sigma],'--r')
grid()
show()


