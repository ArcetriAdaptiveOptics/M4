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
# allinea la parabola , aggiunto in mini_ott/measurements 
meas.ParAlign()

## 
#burst, aggiunto in interf
# ATTENZIONE!! NON ESISTE PIU'!!!!!
# interf.burstAndConvertFrom4DPCTom4OTTpc(1000)

tn=interf.capture(1000)
interf.produce(tn)


##
# fai N misure con delay_run lunga

delay=25
nmeas = 1000
attesa_min=120
attesa_s=attesa_min*60

t1=time.time()
tn = meas.opticalMonitoring(nmeas, delay, attesa_s) 
t2=time.time()-t1

print(tn)
print('measure time: {:.1f} s'.format(t2))

##
# fai N misure con delay_run veloce

meas = Measurements(ott, interf) #crea la classe 
delay=5
nmeas = 100
attesa_s=0

t1=time.time()
tn = meas.opticalMonitoring(nmeas, delay, attesa_s) 
t2=time.time()-t1

print(tn)
print('measure time: {:.1f} s'.format(t2))


fl = th.fileList(tn)
nf = len(fl)
q0=th.averageFrames(0,nmeas-1,fl)
coeff, mat = zernike.zernikeFit(q0, [5,6,7,8])
print(coeff)


##
# faccio misure veloce ogni tot minuti

delay=2
nmeas = 100
attesa_min=60
attesa_s=attesa_min*60
Nset=24

time.sleep(60*60*5)

tn=[]
for ii in range(Nset):
    
    print('start measure {:n}/{:n}'.format(ii+1,Nset))
    tn.append(meas.opticalMonitoring(nmeas, delay, 0))
    print(tn[ii])
    print('wating ...')
    
    time.sleep(attesa_s)

