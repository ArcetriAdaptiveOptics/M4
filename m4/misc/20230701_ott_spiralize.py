import numpy as np
from matplotlib import pyplot as plt
from m4 import main
from astropy.io import fits as pyfits
from m4.configuration import start
from m4.devices.i4d import I4D
from m4.configuration.ott_parameters import Interferometer
conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
ott, interf, dm = start.create_ott(conf)

def spiral_pos(nstep, step):
    p = np.array([0,0])
    pp = []
    for i in range(nstep):
            direct = (-1)**i
            mov = i+1
            for j in range(mov):
                    p = p+np.array([direct,0])
                    pp.append(p)
            for j in range(mov):
                    p = p+np.array([0,direct])
                    pp.append(p)
    pp = np.array(pp)
    plt.plot(pp[:,0],pp[:,1],'-x')
    return pp



def spiralize(p):
    npos = len(p)
    for i in range(npos):
        p0 = ott.parabola.getPosition()
        p1 = p0 + np.array([0,0,0,p[i,0],p[i,1],0])
        print('New Par command:');print(p1)
        ott.parabola.setPosition(p1)
        r0 = ott.referenceMirror.getPosition()
        r1 = r0 + np.array([0,0,0,-2*p[i,0],-2*p[i,1],0])
        print('New RM command:'); print(r1)
        ott.referenceMirror.setPosition(r1)
        



