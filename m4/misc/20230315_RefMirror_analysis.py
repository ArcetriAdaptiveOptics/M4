
from astropy.io import fits as pyfits
from m4.mini_OTT import timehistory as th
from matplotlib import pyplot as plt
from m4.ground import zernike
import numpy as np
from m4 import noise
from m4.configuration import config_folder_names as fold_name
import os


tn0='20230315_155719'

tn1='0002_0000'

fl = th.fileList(os.path.join(tn0,tn1))
nf = len(fl)
q0=th.averageFrames(0,nf-1,fl)