import numpy as np
from matplotlib import pyplot as plt
from m4.configuration import start
from m4.mini_OTT import timehistory as th
from m4.ground import zernike as zern
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
# allinea la parabola , 

tnc = '20230113_102942'

zern2corrf = np.array([0,1,2]) #TipTilt focus
dofidf = np.array([0,1,2])# parpist, ParTip, ParTilt
par0=ott.parabola.getPosition()

tna = main.align_PARAndRM(ott, interf, tnc, zern2corrf, dofidf,n_frames=5)

meas.ParAlign()



##
# fai N misure con delay_run lunga


delay=25
nmeas = 1000
print('aspetto')
attesa_min=120
attesa_s=attesa_min*60
time.sleep(5400) #aspetta tot secondi prima di iniziare
print('parto con le misure')

t1=time.time()
tn = meas.opticalMonitoring(nmeas, delay) 
t2=time.time()-t1

print(tn)
print('measure time: {:.1f} s'.format(t2))

##
# fai N misure con delay_run veloce

meas = Measurements(ott, interf) #crea la classe 
delay=.8
nmeas = 100
print('aspetto')
#time.sleep(5400) #aspetta tot secondi prima di iniziare
print('parto con le misure')

t1=time.time()
tn = meas.opticalMonitoring(nmeas, delay) 
t2=time.time()-t1

print(tn)
print('measure time: {:.1f} s'.format(t2))



## 
#burst

dd = 'D:/M4/Capture/'
tn = ts.now()
print(tn)
interf.burstFramesToSpecificDirectory(dd+tn+'/', 1000)

# convert the frames 
d1 ='D:/M4/Produced/'
interf.convertRawFrames(d1+tn,dd+tn)

# comando da eseguire nella shell con #####=tn
# rsync -av '/home/m4/4d/M4/Produced/#####' '/mnt/m4storage/Data/M4Data/OPTData/OPDImages/'
##
# tn='20230210_123559'
# tau_vector=np.arange(1,200,2)
# fl=th.fileList(tn)
# sf=th.strfun_fl(fl, tau_vector, [1,2,3,4,7,8])
