import numpy as np
from matplotlib import pyplot as plt
from m4.configuration import start

conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
ott, interf, dm = start.create_ott(conf)

from m4 import noise

## burst
tn=[]

tn.append('20230220_163535') # fanBig v3 45°, fanSmall v2 ?° 
tn.append('20230221_163048') # fanBig v3 45°, fanSmall v2 0° 
# tn.append('20230220_170323') # fanBig v3 45°, fanSmall v2 0° 
# tn.append('20230221_184553') # fanBig v3 45°, fanSmall v2 0-5?°
# tn.append('20230221_185331') # fanBig v3 45°, fanSmall v2 0-5?°
# tn.append('20230221_190430') # fanBig v3 45°, fanSmall v2 0-5?°
tn.append('20230222_091309') # fanBig v3 45°, fanSmall v2 0-5?° 

#
tn=[]
tn.append('20230224_112849') # fanGig v1 45°, fanSmall NE v1 0-5?
tn.append('20230224_114327') # fanGig v1 45°, fanSmall NE v2 0-5?
tn.append('20230224_115849') # fanGig v1 45°, fanSmall NE v1 sulla colonna gialla
tn.append('20230224_121353') # fanGig v1 45°, fanSmall E v2 0-5?

# new pad
tn=[]
tn.append('') # fanGig v , fanSmall E v 
tn.append()
tn.append()
tn.append()
tn.append()




## structure funcion on burst
tau_vector = np.arange(1,100,2)
path_series = '/mnt/m4storage/Data/M4Data/OPTData/OPDImages/'

plt.close('all')

for jj in tn:
    dfpath=path_series+jj
    noise.convection_noise(dfpath, tau_vector)

