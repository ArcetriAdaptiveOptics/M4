# script utili per fare analisi su set di misure a 1 min l'una dall'altra per 10 ore
# obiettivo: capire quali condizioni ci permettono di avere una misura "vera" della forma della parabola

from astropy.io import fits as pyfits
from m4.mini_OTT import timehistory as th
from matplotlib import pyplot as plt
from m4.ground import zernike
import numpy as np
from m4 import noise
from m4.configuration import config_folder_names as fold_name
import os

# definisco la maschera ma
nx=425
ny=425
x = np.linspace(1, nx, nx)
y = np.linspace(1, ny, ny)
xv, yv = np.meshgrid(x, y)
x0=215.5;y0=215.5;R=26.5+8;
mask1=(xv-x0)**2+(yv-y0)**2<R**2
x0=215.5;y0=214.5;R=196.5-10;
mask2=(xv-x0)**2+(yv-y0)**2>R**2
ma=mask1 | mask2

## run notturni e diurni presi a condizioni diverse
#60 sec, 1000 misure

tn=[]
tn.append('20230127_163227')
tn.append('20230129_100739')  # day!!
tn.append('20230131_181616')  # cover
tn.append('20230201_183722')  # cover
tn.append('20230202_175230')  # cover, fan
tn.append('20230203_182356')  # cover, fan, box on
tn.append('20230204_160455')  # cover, fan, box on

## run notturni presi alle stesse condizioni
#60 sec, 1000 misure

tn=[]
tn.append('20230203_182356')  # cover, fan, box on
tn.append('20230204_160455')  # cover, fan, box on
tn.append('20230210_181437')  # cover, fan, box on
tn.append('20230213_181037')  # cover, fan, box on

## run diurni 
#60 sec, 1000 misure

tn=[]
tn.append('20230129_100739')  # day!!
tn.append('20230211_120157')  # day, cover, fan, box on
tn.append('20230212_104309')  # day, cover, fan, box on

## run con ventoloni nuovi 
# 2 sec 200 misure
tn=[]
# fanSmall puntata in cima alla torre
tn.append('20230216_143942') # fanBig v3, lontano 
tn.append('20230216_160931') # fanBig v1, lontano
tn.append('20230216_162152') # fanBig v1, medio
tn.append('20230216_163221') # fanBig v1, vicino
tn.append('20230217_091055') # fanBig v3, lontano, 45°, vedo calore scatola 
# fanSmall spostata davanti a fanLarge puntata dritta sullo specchio
tn.append('20230217_103629') # fanBig v3, fanSmall v1
# 1 sec 100 misure
tn=[]
tn.append('20230217_111359') # fanBig v3 45°, fanSmall v0 
tn.append('20230217_111859') # fanBig v3 45°, fanSmall v2 0°
tn.append('20230217_112502') # fanBig v3 45°, fanSmall v2 0°
tn.append('20230217_111002') # fanBig v3 45°, fanSmall v2 30°
tn.append('20230217_112925') # fanBig v3 45°, fanSmall v2 30°

# 0.8 s 500 misure
tn=[]
tn.append('20230217_155624') # fanBig v3 45°, fanSmall v2 ?° 
# 0.8 s 100 misure
tn.append('20230220_170929') # fanBig v3 45°, fanSmall v2 0° 

# 25-30 s 1000 misure 
tn=[]
tn.append('20230218_100353') # fanBig v3 45°, fanSmall v2 ?°  Day
tn.append('20230219_101330') # fanBig v3 45°, fanSmall v2 ?°  Day
tn.append('20230217_175840') # fanBig v3 45°, fanSmall v2 ?°  Night
tn.append('20230219_221356') # fanBig v3 45°, fanSmall v2 ?°  Night
tn.append('20230220_200134') # fanBig v3 45°, fanSmall v2 0°  Night
tn.append('20230221_211019') # fanBig v3 45°, fanSmall v2 0-5° Night 
tn.append('20230223_224400') # fanBig v3 45°, fanSmall v2 0-5° Night 

## 1s 100 misure ogni 30 minuti 
tn=[]
tn.append('20230222_121504')
tn.append('20230222_124904')
tn.append('20230222_132302')
tn.append('20230222_135701')
tn.append('20230222_143059')
tn.append('20230222_150459')
tn.append('20230222_153857')

tn=[]
tn.append('20230224_142807')
tn.append( '20230224_151704')
tn.append( '20230224_160602')
tn.append( '20230224_165500')
tn.append( '20230224_174358')
tn.append( '20230224_183256')
tn.append( '20230224_192154')
tn.append( '20230224_201051')
tn.append( '20230224_205948')
tn.append( '20230224_214845')
tn.append( '20230224_223743')
tn.append( '20230224_232641')
tn.append( '20230225_001538')
tn.append( '20230225_010435')
tn.append( '20230225_015332')
## 2s 100 misure, fanBig v3 45°, fanSmall diverse posizioni 0-5°  
tn=[]
tn.append('20230224_113450') # fanBig v1 45°, fanSmall v1 0-5° NE - pavimento
tn.append('20230224_114947') # fanBig v1 45°, fanSmall v2 0-5° NE - pavimento
tn.append('20230224_120430') # fanBig v1 45°, fanSmall v1 0-5° NE - colonna gialla
tn.append('20230224_140647') # fanBig v1 45°, fanSmall v2 0-5° E - posizione di prima


## cambio pad
tn=[]
tn.append('20230321_210648') # 500 misure, 10s delay, fans ON pos?? 
tn.append('20230321_225557') # 500 misure, 10s delay, fans ON pos?? 
tn.append('20230322_072934') # 500 misure, 10s delay, fans ON pos??
tn.append('20230322_133607') # 100 misure, 2s delay, fans ON pos?? 


## run vari a confronto, copio quì sotto quello che mi interessa confrontare
tn=[] 
tn.append('20230224_140647') # fanBig v1 45°, fanSmall v2 0-5° E - posizione di prima
tn.append('20230322_133607') # 100 misure, 5s delay, fans?? 

##
#per ogni tn calola la media, e fa il plot togliendo zernike di allineamento, 
#poi stampo gli astigmatismi
tlen=len(tn)
plt.close('all')
# idx=[100, 900]
idx=[0, 99]

qm=[]
zernCoef=np.zeros([tlen,2]);
for j in range(tlen):

    np.disp(j)

    fl = th.fileList(tn[j])
    nf = len(fl)

    q0=th.averageFrames(idx[0],idx[1],fl)
    q0.mask=ma
    plt.figure(j)
    coeff, mat = zernike.zernikeFit(q0, [5,6])
    zernCoef[j,:]=coeff
    tmp=th.removeZernike(q0,[1,2,3,4,7,8])
    qm.append(tmp)
    plt.imshow(tmp)
    plt.title(tn[j]+' RMS={:.2e}'.format(np.std(tmp))+'\n z5={:.2e}'.format(zernCoef[j,0])+', z6={:.2e}'.format(zernCoef[j,1])); 
    
plt.colorbar()    

plt.show()

for j in range(tlen):
   print(tn[j]+':  z5={:.2e}'.format(zernCoef[j,0])+', z6={:.2e}'.format(zernCoef[j,1]))
   # print(zernCoef[j,:]/zernCoef[j,0])

##
#per ogni tn calola la media, poi sottrae la media a tutte le singole misure, toglie zernike di allineamento e plotta l'RMS (con runningmean)

plt.close('all')
tlen=len(tn)
idx=[0, 99]

rrr=[]

for j in range(tlen):

    np.disp(j)

    fl = th.fileList(tn[j])
    nf = len(fl)
    q0=th.averageFrames(idx[0],idx[1],fl)

    rr = []
    for i in range(idx[0],idx[1]):
        qi = th.removeZernike(th.frame(i,fl)-q0,[1,2,3,4,7,8])
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
    tmp=th.removeZernike(q0-q1,[1,2,3,4,7,8])
    plt.imshow(tmp); plt.title(tn[j]+' RMS={:.2e}'.format(np.std(tmp))); plt.colorbar()
plt.show()


## 
# ciclo su tutti i tn
# divido ogni set in 10 subset, calcolo  la media all'interno del subset, poi faccio delle differenze tra i vari subset
tlen=len(tn)

plt.close('all')
for j in range(tlen):

    np.disp(j)
    
    fl = th.fileList(tn[j])
    nf = len(fl)
    
    q=[]
    nn=100
    for i in range(10):
        qq=th.removeZernike(th.averageFrames(i*100,i*100+nn-1,fl),[1,2,3,4,7,8])
        q.append(qq)
       # plt.subplot(3,3,i+1)
       # plt.imshow(qq); plt.colorbar()
        
    plt.figure(j+1,figsize=(12,8))
      
    l1=1;l2=2; plt.subplot(2,3,1); 
    mm=np.std(q[l1]-q[l2])
    plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
    l1=1;l2=4; plt.subplot(2,3,2)
    mm=np.std(q[l1]-q[l2])
    plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
    l1=1;l2=7; plt.subplot(2,3,3)
    mm=np.std(q[l1]-q[l2])
    plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
    l1=3;l2=4; plt.subplot(2,3,4)
    mm=np.std(q[l1]-q[l2])
    plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
    l1=3;l2=6; plt.subplot(2,3,5)
    mm=np.std(q[l1]-q[l2])
    plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))
    l1=3;l2=9; plt.subplot(2,3,6)
    mm=np.std(q[l1]-q[l2])
    plt.imshow(q[l1]-q[l2]);plt.colorbar();plt.title(str(l1)+'-'+str(l2)+', RMS={:.2e}'.format(mm))

    plt.suptitle(tn[j]) 

plt.show()

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
        qq=th.removeZernike(th.averageFrames(i*300,i*300+nn-1,fl),[1,2,3,4,7,8])
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

## Calcolo le medie all'interno dei subset e faccio la differenza con tutti gli altri, plottandone anche l'RMS
tlen=len(tn)
sets=range(tlen)
idx=[0, 99]

plt.close('all')

pp=np.ma.zeros([len(sets),len(sets),425,425])
kk=-1
for ii in sets:
    kk=kk+1
    ll=-1
    for jj in sets:
        print(ii); print(jj)        
        ll=ll+1
        fl1= th.fileList(tn[ii]); 
        fl2= th.fileList(tn[jj]); 
        nf1 = len(fl1); nf2=len(fl2)
        q1=th.averageFrames(idx[0],idx[1],fl1); q2=th.averageFrames(idx[0],idx[1],fl2)
        q1.mask=ma; q2.mask=ma
        q=th.removeZernike(q1-q2,[1,2,3,4,7,8])
        pp[kk,ll,:,:]=q   

mis=pp.min()
mas=pp.max()

rr=np.zeros([len(sets),len(sets)])

for ii in range(len(sets)):
    for jj in range(len(sets)):
        rr[ii,jj]=pp[ii,jj,:,:].std()
plt.figure()
plt.imshow(rr); plt.colorbar()


f=plt.figure(figsize=(12,8))
ll=0
for ii in range(len(sets)):
    for jj in range(len(sets)):
        ll=ll+1
        if jj>ii:
            plt.subplot(len(sets),len(sets),ll)
            plt.imshow(pp[ii,jj,:,:])
            # plt.clim(mis/2,mas/2)   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            mm=np.std(pp[ii,jj,:,:]); plt.title('RMS={:.2e}'.format(mm))
            plt.axis('off')
            plt.colorbar()
  #
plt.tight_layout()
plt.figure()
plt.imshow(pp[0,1,:,:])
#plt.clim(mis/2,mas/2)
plt.colorbar()

plt.show()

## visualizza le differenze in sequenza
plt.figure(figsize=(12,8))
for ii in range(len(sets)):
    for jj in range(len(sets)):
        ll=ll+1
        if jj>ii:
            plt.clf()
            plt.imshow(pp[ii,jj,:,:])
            mm=np.std(pp[ii,jj,:,:]); plt.title('RMS={:.2e}'.format(mm))
            plt.axis('off')
            plt.colorbar()
            plt.pause(0.5)
            plt.show()
