from m4.configuration import start
from m4.mini_OTT import timehistory as th
from m4.ground import zernike
from m4.ground import timestamp
import time as tt
from m4 import noise
import os
import numpy as np
from matplotlib import pyplot as plt
from m4.ground import geo
from importlib import reload
from m4.misc import image_registration_lib as imgreg
from astropy.io import fits as pyfits
import glob
from tabulate import tabulate

def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

runa =0
##  define the tracking numbers 
if runa :
    tnlist = ['20231024_155008', '20231024_155130', '20231024_155250', '20231024_155409', '20231024_155532', '20231024_155654', '20231024_155824', '20231024_155947', '20231024_160111', '20231024_160234', '20231024_160356', '20231024_160518', '20231024_160641', '20231024_160806', '20231024_160927', '20231024_161050', '20231024_161214', '20231024_161337', '20231024_161501', '20231024_161625']
    h0 = (130-600)/1000
    h1 = (890-600)/1000
    delta_move = (h1-h0)/len(tnlist)
    
else:
    tn0='20231010_130115'
    tn1 ='20231010_190535'
    tnlist = th.tnscan(tn0, tn1)
    delta_move=0.05
    h0 = 0
    h1 = 295/1000
    



ntn = len(tnlist)

## open average if already produced
qc = []
c_ave = []
c_std = []
for tn_curr in tnlist:
    coeffs=[]
    path = th.foldname.OPD_SERIES_ROOT_FOLDER+'/'
    flist = glob.glob(path+tn_curr+'/20*.fits')
    zlist = [1,2,3,4,5,6,7,8,9, 10, 11]
    for fname in flist:
        tmp_img = th.read_phasemap(fname)
        if len(coeffs) == 0:
            cc,mat=zernike.zernikeFit(tmp_img,zlist)
            mm = tmp_img.mask ==0
            umat = np.linalg.pinv(mat)
        else:
            tmp_data = tmp_img.data[mm]
            tmp_data[np.isnan(tmp_data)]=0
            cc = umat @ tmp_data
        coeffs.append(cc)
    coeffs = np.array(coeffs)
    c_ave.append(coeffs.mean(axis=0)) 
    c_std.append(coeffs.std(axis=0))
    q=imgreg.load_ott(tn_curr, zlist=[1,2,3])
    qc.append(q) 
    
c_ave = np.array(c_ave)
c_std = np.array(c_std)
## Here we compute the zernike coefficient for focus and astigmatism and coma for each sub apertures
#plt.close('all')
plt.figure()
for ii in np.arange(5)+3:
     plt.plot(c_ave[:,ii])
plt.legend(['Z4','Z5','Z6','Z7','Z8']);
plt.title('Global ZModes on Res, from ParCenter upward')

tnpar  = '20231016_124531' #old '20231011_150952'
px_ott = 0.00076
f0 = 0
f1 = 40
par_remapped = imgreg.load_registeredPar(tnpar)
par_filtered = th.comp_filtered_image(par_remapped,  d=px_ott, verbose=True, disp=False, freq2filter=(f0,f1))
par= par_filtered

coeff_ott=[]
coeff_par=[]
coeff_res=[]
k1=0;k2=1
zlist = [1,2,3,4,5,6,7,8,9,10,11]
for i in  np.arange(k1,ntn-k2,1):
    
    #scelgo se usare la parabola rimappata o filtrata
    par_remasked = imgreg.image_remask(par, qc[i])
    res = qc[i]-2*par_remasked
      
    cc,mat=zernike.zernikeFit(qc[i],zlist)
    coeff_ott.append(cc)

    cc,mat=zernike.zernikeFit(2*par,zlist)
    coeff_par.append(cc)
    
    cc,mat=zernike.zernikeFit(res,zlist)
    coeff_res.append(cc)
    
    # faccio il fit degli zernike sulla pupilla grossa invece che su quella piccola. NON FUNZIONA?
    # cc,mat=zernike.zernikeFitAuxmask(res, par_remapped.mask, [1,2,3,4,5,6,7,8])
    # coeff_aux.append(cc)
    

coeff_ott=np.array(coeff_ott)
coeff_par=np.array(coeff_par)
coeff_res=np.array(coeff_res)
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(11,5))
ax1.plot(coeff_ott[:,3]);
ax1.plot(coeff_ott[:,4]);
ax1.plot(coeff_ott[:,5]);
ax1.plot(coeff_ott[:,6]);
ax1.plot(coeff_ott[:,7])#; ax1.set_ylim([-5e-8,5e-8])

ax2.plot(coeff_par[:,3]);
ax2.plot(coeff_par[:,4]);
ax2.plot(coeff_par[:,5]);
ax2.plot(coeff_par[:,6]);
ax2.plot(coeff_par[:,7])#; ax2.set_ylim([-5e-8,5e-8]) 


ax1.set_title('ott'); ax2.set_title('2*par'); 
ax1.legend(['z4','z5','z6','z7','z8']);ax2.legend(['z4','z5','z6','z7','z8']);
plt.show()





num_meas = np.arange(k1, ntn-k2,1)
deltax = num_meas*delta_move+h0
plt.figure()
plt.errorbar(deltax,coeff_res[:,3],yerr=c_std[num_meas,3]/np.sqrt(len(coeff_res)));
plt.errorbar(deltax,coeff_res[:,4],yerr=c_std[num_meas,4]/np.sqrt(len(coeff_res)));
plt.errorbar(deltax,coeff_res[:,5],yerr=c_std[num_meas,5]/np.sqrt(len(coeff_res)));
plt.errorbar(deltax,coeff_res[:,6],yerr=c_std[num_meas,6]/np.sqrt(len(coeff_res)));
plt.errorbar(deltax,coeff_res[:,7],yerr=c_std[num_meas,7]/np.sqrt(len(coeff_res)));
plt.title("Residue Zernike decomposition")
plt.legend(['z4','z5','z6','z7','z8'])

#subaperture ratio
Rpar=0.72
Rmir=0.3
ratio=Rmir/Rpar

a7 = np.mean( coeff_res[:,6])/(ratio**3)
a8 = np.mean( coeff_res[:,7])/(ratio**3)

#coma contribution in astigmatism coefficient
cz4 = np.polyfit(deltax,coeff_res[:,3], 2)
cz5 = np.polyfit(deltax,coeff_res[:,4], 2)
cz6 = np.polyfit(deltax,coeff_res[:,5], 2)
cz7 = np.polyfit(deltax,coeff_res[:,6], 2)
cz8 = np.polyfit(deltax,coeff_res[:,7], 2)
cz9 = np.polyfit(deltax,coeff_res[:,8], 2)
cz10 = np.polyfit(deltax,coeff_res[:,9], 2)
cz11 = np.polyfit(deltax,coeff_res[:,10], 2)
#estimating alpha (from astigmatism)
#alpha = np.atan(np.mean(ca1/ca2* coeff_res[:,5]/coeff_res[:,4] ))


plt.figure()
plt.errorbar(deltax,coeff_res[:,3],yerr=c_std[num_meas,3]/np.sqrt(len(coeff_res)));
plt.errorbar(deltax,coeff_res[:,4],yerr=c_std[num_meas,4]/np.sqrt(len(coeff_res)));
plt.errorbar(deltax,coeff_res[:,5],yerr=c_std[num_meas,5]/np.sqrt(len(coeff_res)));
plt.errorbar(deltax,coeff_res[:,6],yerr=c_std[num_meas,6]/np.sqrt(len(coeff_res)));
plt.errorbar(deltax,coeff_res[:,7],yerr=c_std[num_meas,7]/np.sqrt(len(coeff_res)));
plt.errorbar(deltax,coeff_res[:,8],yerr=c_std[num_meas,8]/np.sqrt(len(coeff_res)));
plt.errorbar(deltax,coeff_res[:,9],yerr=c_std[num_meas,9]/np.sqrt(len(coeff_res)));
plt.errorbar(deltax,coeff_res[:,10],yerr=c_std[num_meas,10]/np.sqrt(len(coeff_res)));
plt.legend(['z4','z5','z6','z7','z8', 'z9', 'z10', 'z11'])
plt.plot(deltax, cz4[2]+cz4[1]*deltax +cz4[0]*deltax**2, 'k--')
plt.plot(deltax, cz5[2]+cz5[1]*deltax +cz5[0]*deltax**2, 'k--')
plt.plot(deltax, cz6[2]+cz6[1]*deltax +cz6[0]*deltax**2, 'k--')
plt.plot(deltax, cz7[2]+cz7[1]*deltax +cz7[0]*deltax**2, 'k--')
plt.plot(deltax, cz8[2]+cz8[1]*deltax +cz8[0]*deltax**2, 'k--')
plt.plot(deltax, cz9[2]+cz9[1]*deltax +cz9[0]*deltax**2, 'k--')
plt.plot(deltax, cz10[2]+cz10[1]*deltax +cz10[0]*deltax**2, 'k--')
plt.plot(deltax, cz11[2]+cz11[1]*deltax +cz11[0]*deltax**2, 'k--')

from arte.utils.zernike_projection_on_subaperture import ZernikeProjectionOnSubaperture
Zm = np.zeros((len(deltax)*(len(zlist)-3), len(zlist)-3))
for ii in np.arange(len(deltax)):
    zpj= ZernikeProjectionOnSubaperture(Rpar,Rmir,deltax[ii],0)
    Zm_tmp = zpj.getProjectionMatrix() #full Negro matrix
    cr=int((len(zlist)-3)*ii)
    cri=int((len(zlist)-3)*(ii+1))
    Zm[cr:cri, :]=np.transpose(Zm_tmp[2:,2:])


u,w,vt = np.linalg.svd(Zm, full_matrices=False) 
ncut = 0
w_inv = 1/w
if ncut >0:
    w_inv[-ncut:]=0
Zm_pinv = vt.T @ np.diag(w_inv) @ u.T
cc=coeff_res[:,3:].flatten()
align_error_coeffs = Zm_pinv @ cc



for sel in np.arange(8)+3:
#sel=3
    plt.figure()
#    plt.plot(deltax,coeff_res[:,sel])
    data2disp = (Zm @ align_error_coeffs).reshape(len(coeff_res), len(zlist)-3)
    plt.plot(deltax,data2disp[:,sel-3]-coeff_res[:,sel], 'k--') 


#aggiungo l'offset al fit
Zm = np.zeros((len(deltax)*(len(zlist)-3), 2*(len(zlist)-3)))
for ii in np.arange(len(deltax)):
    zpj= ZernikeProjectionOnSubaperture(Rpar,Rmir,deltax[ii],0)
    Zm_tmp = zpj.getProjectionMatrix() #full Negro matrix
    cr=int((len(zlist)-3)*ii)
    cri=int((len(zlist)-3)*(ii+1))
    Zm[cr:cri, 0:8]=np.transpose(Zm_tmp[2:,2:])
    Zm[cr:cri, 8:]=np.diag(np.ones(8))


u,w,vt = np.linalg.svd(Zm, full_matrices=False)
cutm = w<1e-6*np.max(w)
w_inv = 1/w
w_inv[cutm]=0
Zm_pinv = vt.T @ np.diag(w_inv) @ u.T
cc=coeff_res[:,3:].flatten()
ccerr=(c_std[:-1,3:]/np.sqrt(len(coeff_res))).flatten()
fit_result = Zm_pinv @ cc
fit_result2 = Zm_pinv @ ccerr
align_error_coeffs = fit_result[0:8]
align_error_coeff_err = np.abs(np.abs(fit_result2[0:8])-np.abs(fit_result[0:8]))
c0 = fit_result[8:]
c0_err = np.abs(np.abs(fit_result2[8:])-np.abs(fit_result[8:]))



plt.figure()
for sel in np.arange(5)+3:
    plt.plot(deltax,coeff_res[:,sel])

for sel in np.arange(5)+3:
    data2disp = (Zm @ fit_result).reshape(len(coeff_res), len(zlist)-3)
    plt.plot(deltax,data2disp[:,sel-3], 'k--')

plt.legend(['z4','z5','z6','z7','z8', 'Negro84'])
plt.title('Zernike drift (Par removed) in RM movement span [0-300 mm]')

tabella=[]
tabella.append(['[nm]','z4','z5','z6','z7','z8', 'z9', 'z10', 'z11'])
print("Coefficients on Large Pupil (from Negro84) [nm]:")
strout=['Negro84']
for number in align_error_coeffs:
    ss = '{:10.2f}'.format(number*1e9) 
    print(ss)
    strout.append(ss)
tabella.append(strout)
print("Coefficients offset (from local bumps) [nm]:")
strout=['Offset']
for number in c0:
    ss = '{:10.2f}'.format(number*1e9)
    print(ss)
    strout.append(ss)
tabella.append(strout)


print(tabulate(tabella))

#errore in %
print(np.abs(fit_result2/fit_result)*100)

