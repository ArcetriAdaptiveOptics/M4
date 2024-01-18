import numpy as np
from matplotlib import pyplot as plt
from m4.configuration import start
from m4.mini_OTT import timehistory as th
from m4.ground import zernike
from m4.misc import par_meas as pp
from m4.ground import timestamp
from astropy.io import fits as pyfits
import time
from m4 import main
from m4 import noise
from m4.mini_OTT.measurements import Measurements
from m4.configuration.ott_parameters import Interferometer
from m4.devices import i4d
import os

conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
ott, interf, dm = start.create_ott(conf)
ts = timestamp.Timestamp()
meas = Measurements(ott, interf) #crea la classe 

#from oaautils import i4d
i4d = i4d.I4D(Interferometer.i4d_IP, Interferometer.i4d_port)

##


data, height, pixel_size_in_microns, width=i4d.getFringeAmplitudeData()

data2d = np.reshape(data, (width, height))
##
plt.close('all')
plt.figure();plt.clf;plt.imshow(data2d);
plt.clim([20,100])
plt.colorbar();

plt.show()

##
location='/home/m4/Desktop'
file_name='testMarker'
fits_file_name = os.path.join(location, file_name)
pyfits.writeto(fits_file_name, data2d)





