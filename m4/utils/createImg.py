'''
@author: cs
'''

import pyfits
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import ion 


    #cmd=np.random.rand(5352)
def createActMap(cmd):
    fitsFileName='/Users/rm/Desktop/Arcetri/M4/ActuatorCoordinates.fits'
    hduList= pyfits.open(fitsFileName)
    co=np.array(hduList[0].data)
    x=co[0][:] 
    y=co[1][:] 
    
    plt.clf()
    plt.scatter(x, y, c=cmd)
    plt.colorbar()
    ion()  #per continuare a lavorare da terminale
    plt.show()
     
    return 
     
    #plt.plot(x,y, 'ro')   