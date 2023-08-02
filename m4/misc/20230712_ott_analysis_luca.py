# script utili per fare analisi su set di misure a tot min l'una dall'altra per tot ore
# obiettivo: siamo in condizioni di misure "stabili"?

from astropy.io import fits as pyfits
from m4.mini_OTT import timehistory as th
from matplotlib import pyplot as plt
from m4.ground import zernike
import numpy as np
from m4 import noise
from m4.configuration import config_folder_names as fold_name
import os

##
tn=[]
tn.append('20230706_192127')

##

#per ogni tn calola la media, poi sottrae la media a tutte le singole misure, toglie zernike di allineamento e plotta l'RMS (con runningmean)

plt.close('all')
tlen=len(tn)
idx=[200, 799]

rrr=[]

for j in range(tlen):

    np.disp(j)

    fl = th.fileList(tn[j])
    nf = len(fl)
    q0=th.averageFrames(idx[0],idx[1],fl)

    rr = []
    for i in range(idx[0],idx[1]):
        qi = th.removeZernike(th.frame(i,fl)-q0,[1,2,3,4])
        rr.append(qi.std())
    
    rrr.append(rr)
    
plt.figure(10)
    
for j in rrr:
    
    plt.subplot(2,1,1)
    plt.plot(th.runningMean(j,10), '-o'); plt.title('RMS-ref, with runningMean'); plt.legend(tn)

    plt.subplot(2,1,2)
    plt.plot(j, '-o'); plt.title('RMS-ref'); plt.legend(tn)

    plt.show()

##
# quì plotto la differenze tra la media sui frame 100-900 e un sub-set più picoolo
tlen=len(tn)
plt.close('all')
for j in range(tlen):

    fl = th.fileList(tn[j])
    nf = len(fl)

    q0=th.averageFrames(0,950,fl)
    q1=th.averageFrames(850,950,fl)
    plt.figure(j)
    tmp=th.removeZernike(q0-q1,[1,2,3,4])
    plt.imshow(tmp); plt.title(tn[j]+' RMS={:.2e}'.format(np.std(tmp))); plt.colorbar()
plt.show()


## 
# ciclo su tutti i tn
# divido ogni set in 10 subset, calcolo  la media all'interno del subset, poi faccio delle differenze tra i vari subset
# tlen=len(tn)
# 
# plt.close('all')
# for j in range(tlen):
# 
#     np.disp(j)
#     
#     fl = th.fileList(tn[j])
#     nf = len(fl)
#     
#     q=[]
#     nn=100
#     for i in range(10):
#         qq=th.removeZernike(th.averageFrames(i*100,i*100+nn-1,fl),[1,2,3,4,7,8])
#         q.append(qq)
#        # plt.subplot(3,3,i+1)
#        # plt.imshow(qq); plt.colorbar()
#         
#     plt.figure(j+1,figsize=(12,8))
#       
#     l1=1;l2=2; plt.subplot(2,3,1); 
#     mm=np.std(q[l1]-q[l2])
#     plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
#     l1=1;l2=4; plt.subplot(2,3,2)
#     mm=np.std(q[l1]-q[l2])
#     plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
#     l1=1;l2=7; plt.subplot(2,3,3)
#     mm=np.std(q[l1]-q[l2])
#     plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
#     l1=3;l2=4; plt.subplot(2,3,4)
#     mm=np.std(q[l1]-q[l2])
#     plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
#     l1=3;l2=6; plt.subplot(2,3,5)
#     mm=np.std(q[l1]-q[l2])
#     plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
#     l1=3;l2=9; plt.subplot(2,3,6)
#     mm=np.std(q[l1]-q[l2])
#     plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
# 
#     plt.suptitle(tn[j]) 
# 
# plt.show()

##
# ciclo su tutti i tn
# divido ogni set in 3 subset, calcolo  la media all'interno del subset, poi faccio delle differenze tra i vari subset
tlen=len(tn)

plt.close('all')
for j in range(tlen):

    np.disp(j)
    
    fl = th.fileList(tn[j])
    nf = len(fl)
    
    q=[]
    nn=300
    for i in range(3):
        qq=th.removeZernike(th.averageFrames(i*300,i*300+nn-1,fl),[1,2,3,4])
        q.append(qq)
       # plt.subplot(3,3,i+1)
       # plt.imshow(qq); plt.colorbar()
        
    plt.figure(j+1,figsize=(12,8))
      
    l1=0;l2=1; plt.subplot(1,3,1); 
    mm=np.std(q[l1]-q[l2])
    plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
    l1=1;l2=2; plt.subplot(1,3,2)
    mm=np.std(q[l1]-q[l2])
    plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
    l1=0;l2=2; plt.subplot(1,3,3)
    mm=np.std(q[l1]-q[l2])
    plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))


    plt.suptitle(tn[j]) 
plt.show()

    