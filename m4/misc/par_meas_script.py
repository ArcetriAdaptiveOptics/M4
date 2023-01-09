#da far girare prima di incominciare le misure
from m4.misc import par_meas as pp
from m4.ground import read_data as rr
from m4.mini_OTT import timehistory as th 
from m4.mini_OTT.measurements import Measurements 
from m4.configuration import start
conf='/mnt/m4storage/Data/SYSCONFData/m4Config.yaml'
ott, interf, dm = start.create_ott(conf)
meas = Measurements(ott,interf)

#acquisition of time series
NMEAS= 200
DELAY = 5 
tn = meas.opticalMonitoring(NMEAS,DELAY)

#acquisition of single (averaged) frame
img1,tn = pp.acq_par(16,1)
cc=pp.get_zern(tn)

fname='/mnt/m4storage/Data/M4Data/OPTData/PARTest/20221220/20221220_125418.fits'
fname='/mnt/m4storage/Data/M4Data/OPTData/PARTest/20221220/20221220_101045.fits'
sname = '20221220_143412'

fold=th.findTracknum(sname)
flist=th.fileList(sname)
zz =th.zernikePlot(flist[0:10])



img =th.read_phasemap(fname)



#PAR meas Zernike.fits file analysis
from astropy.io import fits as pyfits
#import glob
where = '/mnt/m4storage/Data/M4Data/OPTData/OPD_series/'
fold = '20221220_143412'
fname = 'zernike.fits'
file_name = where + fold + '/' + fname
hdul = pyfits.open(file_name)
data = hdul[0].data
time = data[:,0]
coeff = data[:,1::]
ave = np.average(coeff,0)
sstd = np.std(coeff,0)
figure(figsize=(8,5))
plot(coeff[:,4],'-o', label='ast X')
plot(coeff[:,5],'-o', label='ast Y')
legend(); grid(); ylabel('nm RMS'); xlabel('frame #'); title(fold + ' Astigmatism plot')
figure(figsize=(8,5))
plot(coeff[:,6],'-o', label='coma X')
plot(coeff[:,7],'-o', label='coma Y')
legend(); grid(); ylabel('nm RMS'); xlabel('frame #'); title(fold + ' Coma plot')
print('Average Astigmatism =  %.3e ±  %.3e, %.3e ± %.3e  m' %(ave[4], sstd[4], ave[5], sstd[5]))
print('Average Coma =  %.3e ± %.3e, %.3e ± %.3e m' %(ave[6],sstd[6], ave[7], sstd[7]))




import glob
path = '/mnt/m4storage/Data/M4Data/OPTData/PARTest/'
fold = '20221220'
where = path + fold + '/'
D = sorted(glob.glob(where + '*.fits'))
coeff = []
for ii in range(len(D)):
    img = th.read_phasemap(D[ii])
    cc = pp.get_zern(img)
    coeff.append(cc)
coeff = np.array(coeff)





#D = sorted(glob.glob(where + '20221220_*'))
#coeff_series = []
#for ii in range(len(D)):
#    img = th.read_phasemap(D[ii])
#    cc = pp.get_zern(img)
#    coeff_series.append(cc)








