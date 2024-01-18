import numpy as np
from matplotlib import pyplot as plt
from m4.misc import par_meas as pp
from m4.ground import read_data as rr
from m4.mini_OTT import timehistory as th
from m4.ground import zernike as zern
from m4.mini_OTT.measurements import Measurements
from m4 import main
from astropy.io import fits as pyfits
from m4.configuration import start
from m4.devices.i4d import I4D
from m4.configuration.ott_parameters import Interferometer
from m4 import noise
import time
conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
ott, interf, dm = start.create_ott(conf)
meas = Measurements(ott,interf)
phcamfocus = I4D(Interferometer.i4d_IP, Interferometer.i4d_port)

