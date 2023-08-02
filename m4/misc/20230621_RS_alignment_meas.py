import numpy as np
from matplotlib import pyplot as plt
from m4.configuration import start
from m4.mini_OTT import timehistory as th
from m4.ground import zernike
from m4.misc import par_meas as pp
from m4.ground import timestamp
import time
from m4 import main
from m4 import noise
from m4.mini_OTT.measurements import Measurements

conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
ott, interf, dm = start.create_ott(conf)
ts = timestamp.Timestamp()
meas = Measurements(ott, interf) #crea la classe 

##

tn=interf.capture(10)
print('ora produce')
interf.produce(tn)


##
delay=0
nmeas = 25
attesa_s=0

t1=time.time()
tn = meas.opticalMonitoring(nmeas, delay, attesa_s) 
t2=time.time()-t1

print(tn)
print('measure time: {:.1f} s'.format(t2))