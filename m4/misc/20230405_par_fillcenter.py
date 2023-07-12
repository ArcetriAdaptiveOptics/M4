import numpy as np
from m4.mini_OTT import timehistory as th
from m4.ground import geo
from matplotlib import pyplot as plt

def filterFrames(tn,rad1=0.1,rad2=0.02,thr=1.2):
    fl = th.fileList(tn)
    img1= th.frame(1, fl)
    cir = geo.qpupil(np.invert(img1.mask))
    maskdiam = 1.4
    ps = (cir[2]*2/maskdiam)  #pix/m
    inmrad1 = ps*rad1
    inmrad2 = ps*rad2
    mm = geo.draw_mask(np.zeros(shape(img1)),cir[0],cir[1],inmrad1)
    mm = geo.draw_mask(mm,cir[0],cir[1],inmrad2, out=1)
    mmask = -1*mm+1
    st=[]
    for i in range(100):
        img = th.frame(i,fl)
        img = np.ma.masked_array(img, mmask)
        img = th.removeZernike(img, [1,2,3])
        st.append(img.std())    
    st = np.array(st)
    plt.plot(st)
    p = where(st < np.min(st)*thr)[0]
    print('N. elements found:')
    print(len(p)+'/'+print(len(fl)))

    f = []
    for i in p:
        f.append(fl[i])
    ccube = th.cubeFromList(f)
    aveimg = np.ma.mean(ccube, 0)
    return aveimg


