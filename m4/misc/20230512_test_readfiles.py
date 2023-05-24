import numpy as np
from astropy.io import fits as pyfits
from m4.mini_OTT import timehistory as th
tn='20230508_010159'



def myread(fname):
    hduList = pyfits.open(fname)
    masked_array = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))
    return masked_array

base ='/mnt/data/M4/Data/M4Data/OPTData/OPDSeries/'+tn+'/'
flist = th.fileList(tn)
print(len(flist))
v=[]
for i in range(len(flist)):
    q=myread(flist[i])
    v.append(q)
    print(i)


