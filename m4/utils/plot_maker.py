'''
Autors
  - C. Selmi: written in 2019
'''

from astropy.io import fits as pyfits
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import ion
from m4.configuration.ott_parameters import OttParameters


    #cmd=np.random.rand(5352)
    #a=np.ones(892) b=np.zeros(4460) cmd=np.append(a,b)
def createActMap(cmd):
    '''
    Function to plot values assigned to actuators
    in the case of all m4 or one segment

    Parameters
    ----------
        cmd: numpy array
            vector of command for N_ACTS_TOT or N_ACT_SEG actuators

    Returns
    -------
            graphic plot
    '''
    if cmd.shape[0] == OttParameters.N_ACTS_TOT:
        cmd = cmd
    elif cmd.shape[0] == OttParameters.N_ACT_SEG:
        print('Segment number: ')
        input_number = int(input())
        if input_number < OttParameters.N_SEG:
            segment_index = input_number
        else:
            raise OSError('Segment number %s doesnt exists' % input_number)

        command = np.zeros(OttParameters.N_ACTS_TOT)
        for j in range(cmd.shape[0]):
            act = j + (OttParameters.N_ACT_SEG * segment_index)
            command[act] = cmd[j]
        cmd = command

    else:
        raise OSError('The number of actuators chosen is incorrect')

    fits_file_name = OttParameters.M4COORDINATE_ROOT_FOLDER
    hduList = pyfits.open(fits_file_name)
    co = np.array(hduList[0].data)
    input_number = co[0][:]
    y = co[1][:]

    plt.clf()
    plt.scatter(input_number, y, c=cmd, marker=",") #marker dimensione del pixel
    plt.colorbar()
    ion()  #per continuare a lavorare da terminale
    plt.show()
    return

    #plt.plot(input_number,y, 'ro')
