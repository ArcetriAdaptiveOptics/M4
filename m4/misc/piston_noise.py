#cd 'C:\Users\lucao\git'
from m4.configuration import start
from click.termui import pause
#conf='G:\il mio Drive\Lavoro_brera\M4\LucaConfig.yaml'
#ott, interf, dm = start.create_ott(conf)
import glob
import os
from m4.ground import zernike
from matplotlib import pyplot as plt
from matplotlib import *
import matplotlib
import time
import numpy as np
from m4.utils import osutils as osu

def main():
    a='G:/Il mio Drive/Lavoro_brera/PROGETTI/00_M4/Data'
    tn='20210118_233600_noise'
    imgcube=load_datacube(a,tn)
    cube=remZern_datacube(imgcube[1:500,:,:])
    script_confronto_zone(cube)

    return imgcube,cube
    
def load_datacube(a='G:/Il mio Drive/Lavoro_brera/PROGETTI/M4/00_Data',tn='20210118_233600_noise',Nmax=0):
    fl=fileList(tn,a)
    if Nmax>0 & Nmax<len(fl):
        fl=fl[0:Nmax]
    imgcube = osu.createCube(fl)
    
    return imgcube

def remZern_datacube(imgcube,modes=np.array([1,2,3,4])):
    
    cube=np.zeros(shape(imgcube))
    coeff=np.zeros([len(cube),len(modes)])
    
    for j in range(len(imgcube)):
        cube[j,:,:], coeff[j,:]=removeZernike(imgcube[j,:,:],modes)
    
    figure()
    for j in range(len(modes)):    
        plot(coeff[:,j])   
       
    return cube, coeff

def showMovie(imgcube,diff=0):
    
    for j in range(len(imgcube)):
        plt.pause(0.1)
        plt.show()
        plt.imshow(imgcube[j,:,:])
        

    

def script_singlezone(imgcube,N_diffalg=5):
    
    close('all')
     
    pist=np.ma.mean(imgcube[:,:,:], axis=1)
    pist=np.ma.mean(pist, axis=1)
        
    thr=632e-9/4
    pu = signal_unwrap(pist,thr)
    t=np.arange(len(pist))/28.57
    
    figure(1)
    plot(t,pist-pist[0],'-o')
    plot(t,pu,'-o')
    
    figure(3)
    dd= pist[1:]-pist[0:-1]
    plot(t[0:-1],dd,'-o')
    plot(t[0:-1],np.zeros(len(dd))+thr)
    plot(t[0:-1],np.zeros(len(dd))-thr)
      
    fig=plt.figure(4)
    (ax1, ax2) = fig.subplots(2)
    
    template=np.array([1,3,5,7,11,13,15])
    m,sig=find_best_diffalg(pu,template)
    ax1.plot(template,m,'-o')
    ax2.plot(template,sig,'-o')
       
    figure(5)
    s=diffalg(pu,N_diffalg)
    plot(s,'-o')
    
    return s

    
def script_confronto_zone(imgcube):
    
    close('all')
    #tn='20210118_233600_noise' 
    mx=np.array([245,245,300,200]) 
    my=np.array([190,290,240,240])
    
    
    L=5
    
    plt.imshow(imgcube[1,:,:])
    ax = plt.gca()
    
    fig=plt.figure(5)
    (ax1, ax2) = fig.subplots(2)
    
    for j in range(len(mx)):
        
        disp(j)
          
        rect = matplotlib.patches.Rectangle((mx[j]-L, my[j]-L), 2*L, 2*L, linewidth=1, edgecolor='r', facecolor='none')
        ax.add_patch(rect)
        plt.show()
        
        pist=np.ma.mean(imgcube[:,mx[j]-L:mx[j]+L,my[j]-L:my[j]+L], axis=1)
        pist=np.ma.mean(pist, axis=1)
            
        thr=632e-9/4
        pu = signal_unwrap(pist,thr)
        t=np.arange(len(pist))/28.57
        
        figure(2)
        plot(t,pist-pist[0],'-o')
        
        figure(3)
        plot(t,pu,'-o')
        
        figure(4)
        dd= pist[1:]-pist[0:-1]
        plot(t[0:-1],dd,'-o')
        plot(t[0:-1],np.zeros(len(dd))+thr)
        plot(t[0:-1],np.zeros(len(dd))-thr)
        
        
        template=np.array([1,3,5,7,11,13,15])
        m,sig=find_best_diffalg(pu,template)
        ax1.plot(template,m,'-o')
        ax2.plot(template,sig,'-o')
       
       
        figure(6)
        s=diffalg(pu,N=5)
        plot(s,'-o')
        
       
       
#        figure(5)
#        dd= pist[1:]-pist[0:-1]
#        plot(t[0:-1],np.absolute(dd.data),'-o')
#        plot(t[0:-1],np.zeros(len(dd))+thr)
  
        

    


def signal_unwrap(x, thr=632e-9/4, phase = 632e-9/2):
    v = x-x[0]
    npx = np.size(v)
    for i in np.arange(1,npx):
        dv = v[i]-v[i-1]
        if dv > thr:
            v[i] =v[i]-np.abs(phase * (round(dv/phase)))
        if dv < -thr:
            v[i] =v[i]+np.abs(phase * (round(dv/phase)))
    return v


def signal_unwrap_single(x, thr=632e-9/4, phase = 632e-9/4):
    v=x.copy()
    v = v-v[1:]
    npx = np.size(v)
    for i in np.arange(1,npx):
        #dv = v[i]-v[i-1]
        if v[i] > thr:
            v[i] =v[i]-np.abs(phase * (round(v[i]/phase)))
        if v[i] < -thr:
            v[i] =v[i]+np.abs(phase * (round(v[i]/phase)))
    return v


def diffalg(x,N=3):
    
    dim=int(floor(len(x)/N))
    s=np.zeros(dim)
    
    for j in range(dim):

        for k in range(N-1):
            
            s[j]=s[j]+((x[j*N+k+1]-x[j*N+k])*(-1)**k)/2

    
        s[j]=s[j]/(N-1)
    
    return s

def find_best_diffalg(x,template):
    
    m=np.zeros(len(template))
    sig=np.zeros(len(template))
    
    for i in range(len(template)):
        s=diffalg(x,template[i])
        m[i]=np.ma.mean(s)
        sig[i]=np.ma.std(s)
        
    return m,sig


def fileList(tn, fold=None):
    '''
    Parameters
    ----------
    tn: string
        tracking number where to search for the images file list

    Returns
    -------
    lsdir: the list of image files
    fold1: the path where the file is found
    '''
    if fold is not None:
        name = '*.h5'
        addfold ='/hdf5/'
    else:
        
        fold = findTracknum(tn)
        addfold = '/'
        name = '20*'
        if fold == 'OPDImages':
            addfold = '/hdf5/'
            name = 'img*'

    fold1 = fold+'/'+tn+addfold   #to be re-checked at OTT!! 
    lsdir = sorted(glob.glob(fold1+name), key=lambda x: (os.path.basename(x).split('/')[-1].split('.')[0]))
    #lsdir = lsdir[0]

    return lsdir



def removeZernike(ima, modes=np.array([1,2,3,4])):

        coeff, mat = zernike.zernikeFit(ima, modes)
        surf = zernike.zernikeSurface(ima, coeff, mat)
        new_ima = ima-surf
        return new_ima, coeff


# a='G:/Il mio Drive/Lavoro_brera/PROGETTI/M4/Data'
# #tn='20201215_140200_noise'
# tn='20210118_233600_noise'
# #tn='20210120_072200_noise'
#
#
# fl=fileList(tn,a)
# imgcube = th.cubeFromList(fl)
#
# figure(10)
# imshow(imgcube[1,220:240,180:200])
#
# pist=np.ma.mean(imgcube[:,220:240,180:200], axis=1)
# pist=np.ma.mean(pist, axis=1)
#
#
# thr=632e-9/4
# pu = signal_unwrap(pist,thr)
# t=np.arange(len(pist))/28.57
#
# figure(1)
# plot(t,pist-pist[0],'-o')
# plot(t,pu,'-o')
#
#
# figure(2)
# spec,f=th.spectrum(pu,dt=1/28)
# plot(f,spec)
#
# figure(3)
# dd= pist[1:]-pist[0:-1]
# plot(t[0:-1],dd,'-o')
# plot(t[0:-1],np.zeros(len(dd))+thr)
# plot(t[0:-1],np.zeros(len(dd))-thr)
#
# template=[3,5,7,9,13,17]
# m1,sig1=find_best_diffalg(pist,template)
# m2,sig2=find_best_diffalg(pu,template)
#
# figure(4)
# plot(template,sig1)
# plot(template,sig2)
# figure(5)
# template=[3,5,7,9,13,17]
# plot(template,m1)
# plot(template,m2)


#ps = signal_unwrap_single(pist)
