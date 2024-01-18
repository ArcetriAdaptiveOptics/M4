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
from m4.analyzers import requirement_analyzer as req_anal

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
        temp=geo.draw_mask(np.ones(np.shape(im.data)), mark_pos[0,j], mark_pos[1,j], r=4 , out=0)
        mm=np.logical_and(mm,temp)
    
    return mm

##  define the tracking numbers 
# tn0='20231010_130115'
# tn1 ='20231010_190535'
# 
tn0='20231017_092154'# muovendo la truss
tn1 = '20231017_191653'
tnlist = th.tnscan(tn0, tn1)


ntn = len(tnlist)

# open average if already produced
qc = []
for i in tnlist:
    # q = th.openAverage(i)
    # #plt.figure();plt.imshow(q);plt.clim([-10e-8, 10e-8]);plt.colorbar();plt.show()
    # q = th.removeZernike(q,[1,2,3,4])
    q=imgreg.load_ott(i, zlist=[1,2,3])
    qc.append(q) 


## Here we average all the imaeges of the second half of the dataset, superimposing the center of the external mask, this gives an estimation of the flatMirror
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
#plt.figure(); plt.imshow(par_remapped);plt.plot(ottf[1,:],ottf[0,:],'x') ;  plt.clim([-5e-8,5e-8]); plt.show()
flat=[]
coord=[]
k1=36
k2=2
# uso il frame a 0 per calcolare le coordinate al centro
im=qc[40].copy()
xc0, yc0, r, xx, yy =geo.qpupil(qc[40].mask*-1+1)

m1=mask_linear(im,750,766,1200,1228,10)
m2=mask_linear(im,728,1300 ,1200,850,10)
m3 = geo.draw_mask(np.ones(np.shape(im.data)), xc0, yc0, r=70 , out=1)
m4 = marker_mask(im,ottf) 
mm=np.logical_and(m1,m2)
mm=np.logical_and(mm,m3)
mm=np.logical_and(mm,m4)
im.mask[mm==0]=1
#plt.figure(); plt.imshow(th.removeZernike(im,[1,2,3,4,5,6,7,8,9,10,11])); plt.show()

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
plt.figure();plt.imshow(im_flat);plt.clim([-3e-8, 3e-8]); plt.colorbar();plt.title("RMS: {:.2} m".format(np.std(im_flat)));  plt.show();

##
dove="/home/m4/Desktop/immagini"
name="FlatMirror.fits"
fits_file_name = os.path.join(dove, name)
pyfits.writeto(fits_file_name, im_flat.data)
pyfits.append(fits_file_name, im_flat.mask.astype(int))
##
# image=im_flat
# 
# radius_m=0.04
# pixelscale=1/px_ott
# step=1000 # numero di pixel che salta nel for
# n_patches=None
# 
# a2,b2=req_anal.curv_fit_v2(image,pixelscale); print('raggi in km:',np.round(1e-3/a2),np.round(1e-3/b2))
# 
# [req, list_ima, result_vect]=req_anal.patches_analysis(image,radius_m, pixelscale, step , n_patches)
# 
# mm=result_vect[:,0]<80;
# res=result_vect[mm,0]; plt.figure(); plt.hist(res,10); plt.show()
# mm=result_vect[:,1]<80;
# res=result_vect[mm,1]; plt.figure(); plt.hist(res,10); plt.show()

## Here we average all the images without moving them, this is related with the image of the cavity, parabola+relay
plt.close('all')

mm = geo.draw_mask(np.ones(np.shape(par.data)), np.shape(par)[0]/2, np.shape(par)[1]/2, 800, out=1)
par.mask[mm!=0]=1 
plt.figure();plt.imshow(par); plt.show()

ff=9
qq=qc[17-ff:17+ff]
cav=np.ma.mean(qq,axis=0)

plt.figure();plt.imshow(cav);plt.clim([-3e-8, 3e-8]);
plt.figure();plt.imshow(im_flat);plt.clim([-3e-8, 3e-8]); plt.colorbar();plt.title(str(np.std(im_flat)));

plt.show()

##
im=qc[17]

# ritaglio nella posizione del flat
x0, y0, r, xx, yy =geo.qpupil(im.mask*-1+1)
mm = geo.draw_mask(np.ones(np.shape(im.data)), x0, y0, r-5, out=1)
im.mask[mm!=0]=1 
im=im[int(x0)-440:int(x0)+440,int(y0)-440:int(y0)+440]#=im[int(x0)-440:int(x0)+440,int(y0)-440:int(y0)+440]-im_flat
cav=cav[int(x0)-440:int(x0)+440,int(y0)-440:int(y0)+440]

#plt.figure(); plt.imshow(cav); plt.title("cavity");plt.colorbar(); plt.clim([-50e-8,50e-8]); plt.show()
#plt.figure(); plt.imshow(im); plt.title("image");plt.colorbar(); plt.show()
#plt.figure();plt.imshow(im_flat);plt.clim([-3e-8, 3e-8]); plt.colorbar();plt.title(str(np.std(im_flat)));  plt.show();

l=int(0.05/px_ott) #semilato campionamento 


ind=0
# for cx in np.arange(250,590,85):
#     ixy=ix+1
#     for cy in np.arange(250,590,85):
#         iy=iy+1
q243=[]
for cx in np.arange(310,571,130):
     for cy in np.arange(310,571,130):

        
        cav2=cav[cy-l:cy+l,cx-l:cx+l]; cav2=th.removeZernike(cav2,[1,2,3])
        im2=im[cy-l:cy+l,cx-l:cx+l]; im2=th.removeZernike(im2,[1,2,3])
        im_flat2=im_flat[cy-l:cy+l,cx-l:cx+l]; im_flat2=th.removeZernike(im_flat2,[1,2,3])
        res2=im2-cav2
        res2f=im2-cav2-im_flat2
        
        q243.append(imgreg.quick243(res2f,1/px_ott))
        
        #plt.figure(); plt.imshow(qq); plt.colorbar(); plt.show();
        
        # plt.figure(); plt.imshow(cav2); plt.title("cavity, RMS = {:.2}".format(np.std(cav2)));plt.colorbar(); 
        # plt.figure(); plt.imshow(im2); plt.title("image, RMS = {:.2}".format(np.std(im2)));plt.colorbar(); 
        # plt.figure(); plt.imshow(im_flat2); plt.title("flat, RMS = {:.2}".format(np.std(flat2)));plt.colorbar(); 
        # plt.figure(); plt.imshow(res2); plt.title("im-cav, RMS = {:.2}".format(np.std(res2)));plt.colorbar(); 
        # plt.figure(); plt.imshow(res2f); plt.title("im-cav-flat, RMS = {:.2}".format(np.std(res2f)));plt.colorbar(); 
        # plt.show()
        
        fig, axs = plt.subplots(2,3,figsize=(15, 7))
        ax1=axs[0,0]; ax2=axs[0,1]; ax3=axs[0,2]; ax4=axs[1,0]; ax5=axs[1,1]; ax6=axs[1,2];
        # 
        opts = {'vmin': -5e-8, 'vmax': 5e-8, 'edgecolors': 'none'}
        pc1=ax1.pcolormesh(cav2, **opts); fig.colorbar(pc1); ax1.set_title("cavity (par+relay), RMS = {:.2}".format(np.std(cav2)))
        pc2=ax2.pcolormesh(im2, **opts); fig.colorbar(pc2); ax2.set_title("image on OTT, RMS = {:.2}".format(np.std(im2)))
        opts = {'vmin': -2e-8, 'vmax': 2e-8, 'edgecolors': 'none'}
        pc3=ax3.pcolormesh(im_flat2, **opts); fig.colorbar(pc3); ax3.set_title("RM, RMS = {:.2}".format(np.std(im_flat2)))
        pc4=ax4.pcolormesh(res2, **opts); fig.colorbar(pc4); ax4.set_title("im-cav, RMS = {:.2}".format(np.std(res2)))
        pc5=ax5.pcolormesh(res2f, **opts); fig.colorbar(pc5); ax5.set_title("im-cav-RM, RMS = {:.2}".format(np.std(res2f)))
        pc6=ax6.pcolormesh(q243[ind]*2); fig.colorbar(pc6); ax6.set_title("req 243 (WFE)")
        
        fig.suptitle("@px: "+str(cx) +","+str(cy))
        plt.show()
        ind=ind+1
        plt.savefig(os.path.join("/home/m4/Desktop/immagini","@px: "+str(cx) +","+str(cy)))
 
 ## histogram       
q243=np.array(q243)

lin_q243=np.matrix.flatten(q243)*2
plt.figure();plt.hist(lin_q243,15); plt.title("req 243"); plt.xlabel("WFE [m]")
plt.show() 
plt.savefig(os.path.join("/home/m4/Desktop/immagini","hist 243"))
       
q243_sort=np.sort(lin_q243)
n95=int(len(q243_sort)*95/100)
plt.figure();plt.plot(q243_sort); plt.plot(n95,q243_sort[n95],'o'); plt.title("wavefront error at interact scale"); plt.ylabel("WFE [m]"); plt.xlabel("patch id");plt.show() 
plt.savefig(os.path.join("/home/m4/Desktop/immagini","sorted 243"))
##
plt.close('all')

mm = geo.draw_mask(np.ones(np.shape(par.data)), np.shape(par)[0]/2, np.shape(par)[1]/2, 800, out=1)
par.mask[mm!=0]=1 
plt.figure();plt.imshow(par); plt.show()

qq=qc
cav=np.ma.mean(qq,axis=0)

plt.figure();plt.imshow(cav);plt.clim([-3e-8, 3e-8]);

plt.show()

res2=[]
for jj in np.arange(10,60,1):

    im=qc[jj]
    # ritaglio nella posizione del flat
    x0, y0, r, xx, yy =geo.qpupil(im.mask*-1+1)
    mm = geo.draw_mask(np.ones(np.shape(im.data)), x0, y0, r-5, out=1)
    im.mask[mm!=0]=1 

    cav2=th.removeZernike(cav,[1,2,3])
    im2=th.removeZernike(im,[1,2,3])
    im_flat2=th.removeZernike(im_flat,[1,2,3])

    res=im2-cav2; res[int(x0)-440:int(x0)+440,int(y0)-440:int(y0)+440]=res[int(x0)-440:int(x0)+440,int(y0)-440:int(y0)+440]-im_flat2
    res.mask=im.mask
    res2.append(res)

res2=np.array(res2) 
res3= np.ma.mean(res2,axis=0)
res3.mask=cav.mask
plt.figure();plt.imshow(res3); plt.colorbar(); plt.show();
    
        


