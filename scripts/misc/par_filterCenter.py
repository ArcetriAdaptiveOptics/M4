import numpy as np
from matplotlib.pyplot import *
from m4.mini_OTT import timehistory as th
from m4.ground import geo

#filterFrames(tn,0.1,0.03,1.2)

def getPixelScale(tn):
    fl = th.fileList(tn)
    img1= th.frame(1, fl)
    cir = geo.qpupil(np.invert(img1.mask))
    maskdiam = 1.4
    ps = (cir[2]*2/maskdiam)  #pix/m
    return ps



def filterFrames(tn,rad1=0.1,rad2=0.02,thr=None):
    fl = th.fileList(tn)
    img1= th.frame(1, fl)
    cir = geo.qpupil(np.invert(img1.mask))
   # maskdiam = 1.4
   # ps = (cir[2]*2/maskdiam)  #pix/m
    ps = getPixelScale(tn)
    inmrad1 = ps*rad1
    inmrad2 = ps*rad2
    mm = geo.draw_mask(np.zeros(np.shape(img1)),cir[0],cir[1],inmrad1)
    mm = geo.draw_mask(mm,cir[0],cir[1],inmrad2, out=1)
    mmask = -1*mm+1
#    figure()
#    img1 = np.ma.masked_array(img1.data, mmask)
#    imshow(img1)
    st=[]
    for i in range(len(fl)):
        img=th.frame(i,fl)
        img = np.ma.masked_array(img, mmask)
        img=th.removeZernike(img, [1,2,3])
        st.append(img.std())    
    st = np.array(st)
    print('min value = %.3e m, RMS = %.3e m' %(np.min(st),st.std()))
    figure()
    plot(st,'o')
    if thr is not None:
        thr = thr *np.min(st)
        print('Threshold: %.3e m' %(thr))
    else:
        thr = np.min(st)+st.std()
        print('Threshold: min + 1 sigma = .%.3e m' %(thr))
        print(thr)
    plot([0,len(st)],[thr,thr],'--r',label='threshold = %.3e m' %(thr))
    legend()
    grid()
    title(tn)
    p = np.where(st < thr)[0]
    print('N. elements found: %i' %(len(p)))
    f = []
    for i in p:
        f.append(fl[i])
    ccube = th.cubeFromList(f)
    aveimg = np.ma.mean(ccube,0)
    return aveimg, p

def getFilteredCube(tn,p):
    f = []
    fl = th.fileList(tn)
    for i in p:
        f.append(fl[i])
    ccube = th.cubeFromList(f)
    return ccube

def diffAnalysis(tn,p,bbreak=2,zlist=[1,2,3,4,7,8]):
    ccube = getFilteredCube(tn,p)
    aa = int(len(ccube)/bbreak)
    ave1 = np.ma.mean(ccube[0:aa],0)
    ave2 = np.ma.mean(ccube[-aa:-1],0)
    ave1 = th.removeZernike(ave1,zlist)
    ave2 = th.removeZernike(ave2,zlist)
    dd = ave1 - ave2
    return dd




'''
def zernSub(ccube,zlist=[1,2,3,4,7,8]):
   # Differenze a sottrarre le zernike da singole immagini piuttosto che dalla media
    for ii in range(len(ccube)):
        ccube[ii] = th.removeZernike(ccube[ii], zlist)
    return ccube
'''



