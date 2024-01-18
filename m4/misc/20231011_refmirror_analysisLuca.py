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

def mask_linear(im,x1, y1, x2, y2,d):

    
    # d=10
    x = np.arange(1,np.shape(im)[0]+1,1)
    y = np.arange(1,np.shape(im)[1]+1,1)
    xx, yy = np.meshgrid(x, y)
    m=(y2-y1)/(x2-x1)
    #m=(x2-x1)/(y2-y1)
    q=y1-m*x1
    
    temp1=yy-m*xx-q-d/2>0
    temp2=yy-m*xx-q+d/2<0
    mask=np.array(temp1+temp2)
    #plt.figure(); plt.imshow(mask); plt.show()
    
    return mask

def marker_mask(im,mark_pos):
    
    mm=np.ones(np.shape(im))
    for j in np.arange(1,np.shape(mark_pos)[1],1):
        temp=geo.draw_mask(np.ones(np.shape(im.data)), mark_pos[0,j], mark_pos[1,j], r=4 , out=1)
        mm=np.logical_and(mm,temp)
    
    return mm
##  define the tracking numbers 
#tn0='20231010_130115'
#tn1 ='20231010_190535'

# # 
tn0='20231017_092154'# muovendo la truss
tn1 = '20231017_191653'

# 
# sottraendo parabola e allineando con i modi globali
tnlist = ['20231029_150436' ,'20231029_105812' , '20231029_175926' , '20231028_223658', '20231028_204448' ]
ntn = len(tnlist)


tnlist = th.tnscan(tn0, tn1)
ntn = len(tnlist)

## compute and save the average in each folder
# 
# for i in tnlist:
#     th.saveAverage(i)

## compute the average for each flatMirror positions and remove some zernike
# qc = []
# for i in tnlist:
#     fl = th.fileList(i)
#     q = th.averageFrames(0,99,fl)
#     #q = th.removeZernike(q,[1,2,3,4])
#     qc.append(q) 

## open average if already produced
qc = []
for i in tnlist:
    # q = th.openAverage(i)
    # #plt.figure();plt.imshow(q);plt.clim([-10e-8, 10e-8]);plt.colorbar();plt.show()
    # q = th.removeZernike(q,[1,2,3,4])
    q=imgreg.load_ott(i, zlist=[1,2,3])
    qc.append(q) 

#cc,mat=zernike.zernikeFit(qc[4],np.arange(1,11,1)); print(cc[3])


## Here we average all the images superimposing the center of the external mask, this is related with the image of the flatMirror
plt.close('all')
tnpar  = '20231016_124531' #old '20231011_150952'
px_ott = 0.00076 #mm per px
f0 = 0
f1 = 40
par_remapped = imgreg.load_registeredPar(tnpar)
par_filtered = th.comp_filtered_image(par_remapped,  d=px_ott, verbose=True, disp=False, freq2filter=(f0,f1))

par=par_filtered

tnconf="20231013_230000"
#tnconf="20231011_140000"
cgh_image, ott_image, cghf, ottf = imgreg.init_data(tnconf)
plt.figure(); plt.imshow(par_remapped);plt.plot(ottf[1,:],ottf[0,:],'x') ;  plt.clim([-5e-8,5e-8]); plt.show()

#

# uso il frame a 0 per calcolare le coordinate al centro
idc=40 #40
im=qc[idc].copy()
xc0, yc0, r, xx, yy =geo.qpupil(qc[idc].mask*-1+1)

m1=mask_linear(im,750,766,1200,1228,10)
m2=mask_linear(im,728,1300 ,1200,850,10)
m3 = geo.draw_mask(np.ones(np.shape(im.data)), xc0, yc0, r=70 , out=1)
m4 = marker_mask(im,ottf) 
mm=np.logical_and(m1,m2)
mm=np.logical_and(mm,m3)
mm=np.logical_and(mm,m4)
im.mask[mm==0]=1
#plt.figure(); plt.imshow(th.removeZernike(im,[1,2,3,4,5,6,7,8,9,10,11])); plt.show()

#k1=0;k2=0
flat=[]
coord=[]
k1=36;k2=2
for i in  np.arange(k1,ntn-k2,1):
    par2 = imgreg.image_remask(par, qc[i])
    res = qc[i]-2*par2
    
    x0, y0, r, xx, yy =geo.qpupil(qc[i].mask*-1+1)
    coord.append([x0, y0, r])
    

    res.mask[mm==0]=1
    
    f0=res
    f0_cut=f0[int(x0)-440:int(x0)+440,int(y0)-440:int(y0)+440]
    if np.size(f0_cut)==774400:
        flat.append(f0_cut)

coord=np.array(coord)
x=np.arange(0,np.size(coord[:,0]),1)
z1 = np.polyfit(x, coord[:,0], 1); z2 = np.polyfit(x, coord[:,1], 1); z3 = np.polyfit(x, coord[:,2], 1) 
p = np.poly1d(z1); zz1=p(x)
p = np.poly1d(z2); zz2=p(x)

fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(16,6))
ax1.plot(coord[:,0]-zz1);ax2.plot(coord[:,1]-zz2);ax3.plot(coord[:,2])
ax1.set_title('x-linfit(x)'); ax2.set_title('y-linfit(y)'); ax3.set_title('r')
plt.show()

fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(16,6))
ax1.plot(coord[:,0],'-x');ax2.plot(coord[:,1],'-x');ax3.plot(coord[:,2],'-x')
ax1.set_title('x'); ax2.set_title('y'); ax3.set_title('r')
plt.show()

flat2=np.ma.masked_array(flat)
im_flat=np.ma.mean(flat2,axis=0)
im_flat=th.removeZernike(im_flat,[1,2,3,4,5,6,7,8])
plt.figure();plt.imshow(im_flat);plt.clim([-3e-8, 3e-8]); plt.colorbar();plt.title(str(np.std(im_flat)));  plt.show();

# compute quick 243
qq=imgreg.quick243(im_flat*2,1/px_ott);
plt.figure();plt.imshow(qq);plt.clim([-3e-8, 3e-8]);plt.colorbar();plt.show();plt.title(str(np.std(qq)))
qq2=qq[5:15,5:15]
qv = (np.sort(qq2.flatten()));
val=qv[int(np.sum(np.invert(qq2.mask))*0.95)]
plt.figure();plt.plot(qv); plt.title("req243: "+str(val)); plt.show()

#
opts = {'vmin': -3e-8, 'vmax': 3e-8, 'edgecolors': 'none'}
fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(16,6))
pc1=ax1.pcolormesh(im_flat, **opts);ax1.set_aspect('equal','box')
pc2=ax2.pcolormesh(qq, **opts);ax2.set_aspect('equal','box')
ax3.plot(qv)
ax1.set_title('FlatMirror, RMS: {:.2}'.format(im_flat.std())); ax2.set_title('quick 243'); ax3.set_title('req243: {:.2}'.format(val))
fig.colorbar(pc1); fig.colorbar(pc2)

#fig.suptitle('')
plt.show()
#

## media di tanti ritagli presi lungo il diametro orizzontale dell'immagine del flattone = errore di calibrazione 243

l=int(0.05/px_ott) #semilato campionamento 
step=int(0.01/px_ott) #step di campionamento
start_px=200
stop_px=700

cutx=[]
cuty=[]
cutxy=[]
for i in np.arange(start_px,stop_px,step):
    cutx.append(im_flat[450-l:450+l,i:i+l*2])
    cuty.append(im_flat[i:i+l*2,450-l:450+l])
    cutxy.append(im_flat[450-l:450+l,i:i+l*2])
    cutxy.append(im_flat[i:i+l*2,450-l:450+l])
    
cutx=np.ma.masked_array(cutx)
patchx=np.ma.mean(cutx,axis=0)
cuty=np.ma.masked_array(cuty)
patchy=np.ma.mean(cuty,axis=0)
cutxy=np.ma.masked_array(cutxy)
patchxy=np.ma.mean(cutxy,axis=0)

opts = {'vmin': -0.8e-8, 'vmax': 0.8e-8, 'edgecolors': 'none'}
fig, (ax1, ax2, ax3) = plt.subplots(1,3)
pc1=ax1.pcolormesh(patchx, **opts); fig.colorbar(pc1); ax1.set_aspect('equal','box')
pc2=ax2.pcolormesh(patchy, **opts); fig.colorbar(pc2); ax2.set_aspect('equal','box')
pc3=ax3.pcolormesh(patchxy, **opts); fig.colorbar(pc3); ax3.set_aspect('equal','box')
ax1.set_title('shift along x, RMS: {:.3}'.format(patchx.std()))
ax2.set_title('shift along y, RMS: {:.3}'.format(patchy.std()))
ax3.set_title('shift along x and y, RMS: {:.3}'.format(patchxy.std()))
fig.tight_layout()
plt.show()


## Here we average all the images without moving them, this is related with the image of the parabola. forse
plt.close('all')
# 
tnpar  = '20231016_124531' #old '20231011_150952' #parabola misurata è in 20230428_164257
px_ott = 0.00076
f0 = 0
f1 = 140
par_remapped = imgreg.load_registeredPar(tnpar)
#par_remapped =th.removeZernike(par_remapped ,[1,2,3,4,5,6,7,8])
par_filtered = th.comp_filtered_image(par_remapped,  d=px_ott, verbose=True, disp=False, freq2filter=(f0,f1))
par_dif=par_remapped-par_filtered

par=par_filtered
mm = geo.draw_mask(np.ones(np.shape(par.data)), np.shape(par)[0]/2, np.shape(par)[1]/2, 800, out=1)
par.mask[mm!=0]=1 
plt.figure();plt.imshow(par); plt.show()

#tentativo di correggere l'allineamento della parabola mettendo a zero il centro carotato delle dimensioni del flattone
# cc,mat=zernike.zernikeFit(par,np.arange(1,11,1)); print(cc)
# mm = np.where(par.mask == 0)
# fuoco = np.zeros((par.shape[0], par.shape[1]))
# fuoco[mm] = np.dot(mat[:,3], 20e-9)
# par=par-(np.sign(cc[3])*fuoco)
# 
# par_center=par.copy()
# par_center.mask=qc[0].mask
# cc2,mat2=zernike.zernikeFit(par_center,np.arange(1,11,1))
# # 
# fuoco[mm] = np.dot(mat[:,3], .5e-9)
# while np.abs(cc2[3])>1e-9:
#     par=par-(np.sign(cc2[3])*fuoco)
#     par_center=par.copy()
#     par_center.mask=qc[0].mask
#     cc2,mat2=zernike.zernikeFit(par_center,np.arange(1,11,1))
#     print(cc2[3])
#     
# fuoco[mm] = np.dot(mat[:,3], .1e-9)
# while np.abs(cc2[3])>1e-11:
#     par=par-(np.sign(cc2[3])*fuoco)
#     par_center=par.copy()
#     par_center.mask=qc[0].mask
#     cc2,mat2=zernike.zernikeFit(par_center,np.arange(1,11,1))
#     print(cc2[3])
#     
# fuoco[mm] = np.dot(mat[:,3], .001e-9)
# while np.abs(cc2[3])>5e-13:
#     par=par-(np.sign(cc2[3])*fuoco)
#     par_center=par.copy()
#     par_center.mask=qc[0].mask
#     cc2,mat2=zernike.zernikeFit(par_center,np.arange(1,11,1))
#     print(cc2[3])


# erode some pixel from the external mask to avoid bad reconstruction
temp=[]
RMS=[]
RMS2=[]
jj=0
for i in  np.arange(0,ntn,1):
    
    temp.append(qc[i].copy())
    
    x0, y0, r, xx, yy =geo.qpupil(temp[jj].mask*-1+1)
    mm = geo.draw_mask(np.ones(np.shape(temp[jj].data)), x0, y0, r-5, out=1)
    temp[jj].mask[mm!=0]=1 
    #temp[jj]=th.removeZernike(qc[i],[1,2,3,4,5,6,7,8])
    
    #sottraggo l'immagine del flat calcolata prima
    temp[jj][int(x0)-440:int(x0)+440,int(y0)-440:int(y0)+440]=temp[jj][int(x0)-440:int(x0)+440,int(y0)-440:int(y0)+440]-im_flat
    
    res=temp[jj]-2*par
    res2=th.removeZernike(res,[1,2,3,4,5,6,7,8,9,10])
    RMS.append(res.std())
    RMS2.append(res2.std())
    
    

    jj=jj+1
        
qq=np.ma.masked_array(temp)

plt.figure(); plt.plot(RMS, label=('2*par subtracted'));plt.plot(RMS2, label=('10 Zernike subtracted')); plt.legend();plt.grid(); plt.show(); plt.title("RMS of the residuals moving from the center to the outer regions")



##

#im=np.ma.mean(qq,axis=0)
im = qq[0]

x0, y0, r, xx, yy =geo.qpupil(im.mask*-1+1)
im = im[int(x0)-(int(r)+10):int(x0)+(int(r)+10),int(y0)-(int(r)+10):int(y0)+(int(r)+10)]
par2 = par[int(x0)-(int(r)+10):int(x0)+(int(r)+10),int(y0)-(int(r)+10):int(y0)+(int(r)+10)]
par2.mask=im.mask

# remove zernike
im=th.removeZernike(im,[1,2,3,4,5,6,7,8])
par2=th.removeZernike(par2,[1,2,3,4,5,6,7,8])
#

res=im-2*par2

fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(15, 3.5))

opts = {'vmin': -5e-8, 'vmax': 5e-8, 'edgecolors': 'none'}
pc1=ax1.pcolormesh(im, **opts); fig.colorbar(pc1); ax1.set_title("par from ott images, RMS: {:.3}".format(im.std()))
pc2=ax2.pcolormesh(2*par2, **opts); fig.colorbar(pc2); ax2.set_title("par remapped, RMS: {:.3}".format((2*par2).std()))
pc3=ax3.pcolormesh(res, **opts); fig.colorbar(pc3); ax3.set_title("res, RMS: {:.3}".format(res.std()))
plt.show()

fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(15, 3.5))




## Here we compute the zernike coefficient for focus and astigmatism for each sub apertures
plt.close('all')

# tnpar  = '20231016_124531' #old '20231011_150952'
# px_ott = 0.00076
# f0 = 0
# f1 = 40
# par_remapped = imgreg.load_registeredPar(tnpar)
# par_filtered = th.comp_filtered_image(par_remapped,  d=px_ott, verbose=True, disp=False, freq2filter=(f0,f1))

coeff_ott=[]
coeff_par=[]
coeff_res=[]
k1=40;k2=0
for i in  np.arange(k1,ntn-k2,1):
    
    #scelgo se usare la parabola rimappata o filtrata
    par = imgreg.image_remask(par, qc[i])
    res = qc[i]-2*par
      
    cc,mat=zernike.zernikeFit(qc[i],[1,2,3,4,5,6,7,8])
    coeff_ott.append(cc)

    cc,mat=zernike.zernikeFit(2*par,[1,2,3,4,5,6,7,8])
    coeff_par.append(cc)
    
    cc,mat=zernike.zernikeFit(res,[1,2,3,4,5,6,7,8])
    coeff_res.append(cc)
    
    # faccio il fit degli zernike sulla pupilla grossa invece che su quella piccola. NON FUNZIONA?
    # cc,mat=zernike.zernikeFitAuxmask(res, par_remapped.mask, [1,2,3,4,5,6,7,8])
    # coeff_aux.append(cc)
    

coeff_ott=np.array(coeff_ott)
coeff_par=np.array(coeff_par)
coeff_res=np.array(coeff_res)
fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(16,5))
ax1.plot(coeff_ott[:,3]);ax1.plot(coeff_ott[:,4]);ax1.plot(coeff_ott[:,5]);ax1.plot(coeff_ott[:,6]);ax1.plot(coeff_ott[:,7])#; ax1.set_ylim([-5e-8,5e-8])
ax2.plot(coeff_par[:,3]);ax2.plot(coeff_par[:,4]);ax2.plot(coeff_par[:,5]);ax2.plot(coeff_par[:,6]);ax2.plot(coeff_par[:,7])#; ax2.set_ylim([-5e-8,5e-8]) 
ax3.plot(coeff_res[:,3]);ax3.plot(coeff_res[:,4]);ax3.plot(coeff_res[:,5]);ax3.plot(coeff_res[:,6]);ax3.plot(coeff_res[:,7])#; ax3.set_ylim([-5e-8,5e-8])
ax1.set_title('ott'); ax2.set_title('2*par'); ax3.set_title('res=ott-2*par')
ax1.legend(['z4','z5','z6','z7','z8']);ax2.legend(['z4','z5','z6','z7','z8']);ax3.legend(['z4','z5','z6','z7','z8']); fig.suptitle('Zernike coefficient for different positions in the tower')
plt.show()

# ## Here we compute the zernike coefficient for focus and astigmatism for each sub apertures, più altre cose strane
# #plt.close('all')
# 
# # tnpar  = '20231016_124531' #old '20231011_150952'
# # px_ott = 0.00076
# # f0 = 0
# # f1 = 40
# # par_remapped = imgreg.load_registeredPar(tnpar)
# # par_filtered = th.comp_filtered_image(par_remapped,  d=px_ott, verbose=True, disp=False, freq2filter=(f0,f1))
# 
# cc0,mat0=zernike.zernikeFit(par,np.arange(1,5,1));
# mm0 = np.where(par.mask == 0)
# fuoco = np.zeros((par.shape[0], par.shape[1]))
# 
# coeff_ott=[]
# coeff_par=[]
# coeff_res=[]
# coeff_ott2=[]
# coeff_par2=[]
# coeff_res2=[]
# for i in  np.arange(0,ntn-1,5):
#     
#     par2 = imgreg.image_remask(par, qc[i])
#     res = qc[i]-2*par2
#       
#     cc1,mat1=zernike.zernikeFit(qc[i],[1,2,3,4,5,6])
#     coeff_ott.append(cc1)
# 
#     cc2,mat2=zernike.zernikeFit(2*par2,[1,2,3,4,5,6])
#     coeff_par.append(cc2)
#     
#     cc3,mat3=zernike.zernikeFit(res,[1,2,3,4,5,6])
#     coeff_res.append(cc3)
#     
# 
#     fuoco[mm0] = np.dot(mat0[:,3], cc3[3])*4
#     plt.figure();plt.imshow(fuoco);plt.colorbar();plt.show()
#     fuoco2=imgreg.image_remask(par, qc[i])
#     ott=qc[i]-fuoco2
#     res = ott-2*par2
#       
#     cc4,mat4=zernike.zernikeFit(ott,[1,2,3,4,5,6])
#     coeff_ott2.append(cc4)
# 
#     cc5,mat5=zernike.zernikeFit(2*par2,[1,2,3,4,5,6])
#     coeff_par2.append(cc5)
#     
#     cc6,mat6=zernike.zernikeFit(res,[1,2,3,4,5,6])
#     coeff_res2.append(cc6)
#     
#     # faccio il fit degli zernike sulla pupilla grossa invece che su quella piccola. NON FUNZIONA?
#     # cc,mat=zernike.zernikeFitAuxmask(res, par_remapped.mask, [1,2,3,4,5,6,7,8])
#     # coeff_aux.append(cc)
#     
# 
# coeff_ott=np.array(coeff_ott)
# coeff_par=np.array(coeff_par)
# coeff_res=np.array(coeff_res)
# fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(16,5))
# ax1.plot(coeff_ott[:,3]);ax1.plot(coeff_ott[:,4]);ax1.plot(coeff_ott[:,5]); ax1.set_ylim([-5e-8,5e-8])
# ax2.plot(coeff_par[:,3]);ax2.plot(coeff_par[:,4]);ax2.plot(coeff_par[:,5]); ax2.set_ylim([-5e-8,5e-8]) 
# ax3.plot(coeff_res[:,3]);ax3.plot(coeff_res[:,4]);ax3.plot(coeff_res[:,5]); ax3.set_ylim([-5e-8,5e-8])
# ax1.set_title('ott'); ax2.set_title('2*par'); ax3.set_title('res=ott-2*par')
# ax1.legend(['z4','z5','z6']);ax2.legend(['z4','z5','z6']);ax3.legend(['z4','z5','z6']); fig.suptitle('Zernike coefficient for different positions in the tower')
# plt.show()
# 
# 
# coeff_ott2=np.array(coeff_ott2)
# coeff_par2=np.array(coeff_par2)
# coeff_res2=np.array(coeff_res2)
# fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(16,5))
# ax1.plot(coeff_ott2[:,3]);ax1.plot(coeff_ott2[:,4]);ax1.plot(coeff_ott2[:,5]); ax1.set_ylim([-5e-8,5e-8])
# ax2.plot(coeff_par2[:,3]);ax2.plot(coeff_par2[:,4]);ax2.plot(coeff_par2[:,5]); ax2.set_ylim([-5e-8,5e-8]) 
# ax3.plot(coeff_res2[:,3]);ax3.plot(coeff_res2[:,4]);ax3.plot(coeff_res2[:,5]); ax3.set_ylim([-5e-8,5e-8])
# ax1.set_title('ott'); ax2.set_title('2*par'); ax3.set_title('res=ott-2*par')
# ax1.legend(['z4','z5','z6']);ax2.legend(['z4','z5','z6']);ax3.legend(['z4','z5','z6']); fig.suptitle('Zernike coefficient for different positions in the tower')
# plt.show()