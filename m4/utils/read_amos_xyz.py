import scipy.stats as stats
import scipy.fft
import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np
from arte.utils import surface_fitter as sft

def loadAmos(fxyz):
    wave =None
    wedge = None
    pixsc = None
    count = 0
    with  open(fxyz,'r') as fid:
        while True:
            count +=1
            cline = fid.readline()
            if not cline:
                break
            if (count <15) | ("No Data" in cline) | ("#" in cline):
                #print("Read Line{}: {}".format(count, cline.strip()))
                if count == 4:
                    d1, d2, sx, sy = np.array(cline.split(" "), dtype=int)
                    img = np.zeros([sx,sy])
                    mask = np.zeros([sx,sy])
                if count == 8:
                    dummy, wedge, wave, d1, d2, d3, pixsc, d4 = np.array(cline.split(" "), dtype=float)
                    #FORCED HERE!!!!!
                    wave=1e-6

            else:
                cx,cy,val=map (lambda x:float(x),cline.split(" "))
                cx = int(cx)
                cy = int(cy)
                mask[cx,cy]=1
                img[cx,cy]=val*wave
    sf = sft.SurfaceFitter(np.ma.masked_array(img, mask=(mask==0)))
    sf.fit(np.arange(36)+1)
    zcoeffs = sf.coeffs()
    umat = sf.umat()
    params = sf._coords
    print("pixelscale {} pix/m" .format(pixsc))
    print("wedge {} pix/m" .format(wedge))
    print("wave {} pix/m" .format(wave))
    print("sx {} pix/m" .format(sx))
    print("sy {} pix/m" .format(sy))
    print("params cx,cy,r  {} pix" .format(params))
    print("mean check: mean {}, Z1 {}".format(np.mean(img[mask==1]),zcoeffs[0]))

    return {'img':img, 'mask':mask,'sx':sx,'sy':sy, 'wave':wave, 'params':params,
            'wedge':wedge,'pixsc':pixsc, 'zcoeffs':zcoeffs, 'umat':umat ,'name':fxyz}

