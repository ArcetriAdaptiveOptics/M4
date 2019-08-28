'''
@author: cs
'''

import pyfits
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import ion 

nActsTot= 5352
nActSeg= 892

    #cmd=np.random.rand(5352)
    #a=np.ones(892) b=np.zeros(4460) cmd=np.append(a,b)
def createActMap(cmd):
    if cmd.shape[0]== nActsTot:
        cmd=cmd
    elif cmd.shape[0]== nActSeg:
        print('Segment number: ')
        x = int(input())
        if x < 6:
            segmentIndex = x
        else:
            raise OSError('Segment number %s doesnt exists' % x)
        
        command=np.zeros(nActsTot)
        for j in range(cmd.shape[0]):
            act= j + (nActSeg* segmentIndex)
            command[act]=cmd[j]
        cmd= command
        
    else:
        raise OSError('The number of actuators chosen is incorrect')
    
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