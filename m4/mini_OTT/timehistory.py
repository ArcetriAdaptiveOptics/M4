import os
import glob
import numpy as np
import jdcal
from astropy.io import fits as pyfits
from m4.ground import read_data
from m4.ground import zernike
from m4.ground.read_data import InterferometerConverter
from matplotlib.pyplot import *
import psutil
import scipy.fft
import scipy.stats as stats
ic = InterferometerConverter()
#a= '/mnt/data/M4/Data/M4Data/OPTData/'
a='/mnt/m4storage/Data/M4Data/OPTData/'
tn = '20210425_085242'   #fits
tn  ='20210429_224400_noise'



class TimeHist():
    '''
    Class to 
    HOW TO USE IT::


    '''

    def __init__(self, tn):
        """The constructor """
        self.tracknum = tn

        self._fold  = findTracknum(tn)
        self._path = a+ self._fold
        self._list = fileList(tn)

    def frame(self, id):

        return frame(id, self._list)

    def averageFrames(self, start, stop):
    	
        return averageFrames(start, stop, self._list)



    #@staticmethod
    #def _storageFolder():
    #    """ Creates the path where to save data"""
    #    return fold_name.OPDSERIES



def findTracknum(tn):
    '''
    Parameters
    ----------
    tn: string
        tracking number to be searched in the data folder

    Returns
    -------
    result: the specific data folder where the tracknum is found
    '''

    #a= '/mnt/data/M4/Data/M4Data/OPTData/'
    lsdir = os.listdir(a)
    for i in lsdir:
        b = a+i
        z = os.listdir(b)
        check = False
        for j in z:
            check = (j == tn)
            if check == True:
                result = i
                return result

def _sortFunc4D(elem):
    iid = os.path.basename(elem)[:-3]
    iis=('%5.5i' % int(iid))
    return iis 

def fileList(tn, fold=None, name=None):
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
        if name is None:
            raise
        tn=''
        addfold='/'
        #name = '*.4D'
        #addfold ='/hdf5/'
        fold1 = fold+'/'+tn+addfold
    else:
        
        fold = findTracknum(tn)
        addfold = '/'
        dirs = os.listdir(a+'/'+ fold+'/'+tn)
        if dirs[0] == 'hdf5': 
            addfold = '/hdf5/'
            name = 'img*.h5'
        elif dirs[0][-3:] == '.4D':
            name = '*.4D'
        else:
            name = '20*.fits'
        fold1 =a+'/'+ fold+'/'+tn+addfold   #to be re-checked at OTT!! 
            


    print(fold1+name)
    lsdirs = glob.glob(fold1+name)
    if name == "*.4D":
        lsdirs.sort(key=_sortFunc4D)
        lsdirs.sort(key=_sortFunc4D)
    else:
        lsdirs.sort()
    return lsdirs

def read_phasemap(filename, thefold = None):

    #if thefold is not None:
    #    the = filename.split(thefold)[1]

    thetype = filename.split('.')[1]
    if thetype == 'fits':
        with pyfits.open(filename) as hduList:
            img = hduList[0].data
            mask = hduList[1].data
        #mask = np.zeros(img.shape, dtype=np.bool)
        #mask[np.where(img == img.max())] = True
        img = np.ma.masked_array(img, mask)
        #hduList.close()

    if thetype == '4D':
        img = ic.fromPhaseCam6110(filename)

    if thetype == 'h5':
        img = ic.fromPhaseCam4020(filename)

    return img

def averageFrames(first, last, fileList, thresh=None, fsel=None):
    '''
    Parameters
    ----------
    first: first item 
        tracking number where to search for the images file list

    Returns
    -------
    lsdir: the list of image files
    fold1: the path where the file is found
    '''
    if fsel is None:
        fsel = np.arange(first, last+1)
    
    #imcube = cubeFromList(fileList[x for x in flist])
    imcube = cubeFromList([fileList[x] for x in fsel])
    if thresh is None:
        aveimg = np.ma.mean(imcube, axis=0)
        
    else:
        img = imcube[0].data*0
        mmask = imcube[0].mask
        mysize = imcube[0].compressed()
        nn = 0
        for i in imcube:
            if i.data.compressed.size > 1:
                nn +=1
                img += i.data
                mmask = np.ma.mask_or(i.mask, mmask)
            
        img = img/nn
        image = np.ma.masked_array(img, mask=mmask)
        
    return aveimg

def removeZernike(ima, modes=np.array([1,2,3,4])):

        coeff, mat = zernike.zernikeFit(ima, modes)
        surf = zernike.zernikeSurface(ima, coeff, mat)
        new_ima = ima-surf
        return new_ima
        
def zernikePlot(mylist, modes=np.array(range(1,11))):
    mytype = type(mylist)
    if mytype is list:
        imgcube = cubeFromList(mylist)
        
    if mytype is np.ma.core.MaskedArray:
        imgcube = mylist
        
    zlist = []
    for i in range(len(imgcube)):
        print(i)
        coeff, mat = zernike.zernikeFit(imgcube[i], modes)
        zlist.append(coeff)
        
    zcoeff = np.array(zlist)
    zcoeff = zcoeff.T
    return zcoeff
    
def runningDiff(tn, gap=2):
    llist =fileList(tn)
    nfile = len(llist)
    npoints = int(nfile/gap)-2
    slist=[]
    for i in range(0,npoints):
        #print(i*gap)
        #print(i*gap+1)
        q0 = frame(i*gap,llist)
        q1 = frame(i*gap+1,llist)
        diff = q1-q0
        diff = removeZernike(diff)
        slist.append(diff.std())
    svec = np.array(slist)
    return svec
    
def frame(id, mylist):

    mytype = type(mylist)
    if mytype is list:
        img = read_phasemap(mylist[id])
        
    if mytype is np.ma.core.MaskedArray:
        img = mylist[id]
    
    return img
    
def spectrum(signal, dt=1, show=None):
    # see: https://numpy.org/doc/stable/reference/generated/numpy.angle.html?highlight=numpy%20angle#numpy.angle  for the phase spectrum
    nsig = signal.shape
    if np.size(nsig) ==1:
        thedim = 0

    else:
        thedim = 1
 
    if np.size(nsig) ==1:
        spe  = np.fft.rfft(signal,  norm='ortho')
        nn   = np.sqrt(spe.shape[thedim])   #modRB 


    else:
        spe  = np.fft.rfft(signal, axis=1, norm='ortho')
        nn   = np.sqrt(spe.shape[thedim])   #modRB 
    spe  = (np.abs(spe)) / nn
    freq = np.fft.rfftfreq(signal.shape[thedim], d=dt)
    
    if np.size(nsig) ==1:
        spe[0] = 0
    else:
        spe[:,0] = 0
    if show is not None:
        for i in range(0,nn):
            plot(freq, spe[i,:])
    return spe, freq
        
def cubeFromList(fileList):
    image_list = []
    for i in fileList:
        ima = read_phasemap(i)
        image_list.append(ima)

    image_list = np.ma.masked_array(image_list)
    #print(psutil.Process().open_files())
    return image_list


def timevec(tn):
    fold = findTracknum(tn)
    flist = fileList(tn)
    nfile = len(flist)
    if fold == 'OPDImages':
        tspace = 1./28.57
        timevec = range(nfile)*tspace
    
    
    if fold == 'OPD_series':
        timevec = []
        for i in flist:
            pp = i.split('.')[0]
            tni = pp.split('/')[-1]
            y=tni[0:4]
            mo = tni[4:6]
            d = tni[6:8]
            h = float(tni[9:11])
            mi = float(tni[11:13])
            s = float(tni[13:15])
            jdi=sum(jdcal.gcal2jd(y, mo, d))+h/24+mi/1440+s/86400
            timevec.append(jdi)
        timevec=np.array(timevec)
    
        
    return timevec
        
def track2jd(tni):
    y, mo, d, h, mi, s = track2date(tni)
    jdi=sum(jdcal.gcal2jd(y, mo, d))+h/24+mi/1440+s/86400  
    return jdi
    
def track2date(tni):
    y=tni[0:4]
    mo = tni[4:6]
    d = tni[6:8]
    h = float(tni[9:11])
    mi = float(tni[11:13])
    s = float(tni[13:15])
    return y, mo, d, h, mi, s
    
def runningMean(vec, npoints):
    
    return np.convolve(vec, np.ones(npoints), 'valid') / npoints       
        
        
def strfunct(vect, gapvect):
    '''
    vect shall be npoints x m
    the strfunct is calculate m times over the npoints time series
    returns stf(n_timeseries x ngaps)
    '''
    nn      = np.shape(vect)
    maxgap  = np.max(gapvect)
    ngap    = len(gapvect)
    n2ave   = int(nn/(maxgap))-1 # or -maxgap??
    jump    = maxgap
    st      = np.zeros(ngap)
    for j in range(ngap):
        tx = []
        for k in range(n2ave):
            print('Using positions:')
            print(k*jump,k*jump+gapvect[j])
            tx.append( (vect[k*jump]-vect[k*jump+gapvect[j]])**2)
        st[j]=np.mean(np.sqrt(tx))
    return st

# def strfun_fl(fl, gapvect, ZernModes=None):
#     '''
#     the strfunct is calculate m times over the npoints time series
#     returns stf(n_timeseries x ngaps)
#     '''
#     d=len(fl)
#     
#     cube = []
#     for i in fl:
#         ima = read_phasemap(i)
#         cube.append(ima)
# 
#     cu1 = np.ma.masked_array(cube)
#     
#     cu2=[]
#     if ZernModes is None:
#         cu2=cu1
#     else:
#         for jj in range(d):
#             cu2.append(removeZernike(cu1[jj,:,:],ZernModes))
#     
#     disp('imagecube created and zernike removed')
#     
#     nn      = len(cu2)
#     maxgap  = np.max(gapvect)
#     ngap    = len(gapvect)
#     n2ave   = int(nn/(maxgap))-1 # or -maxgap??
#     jump    = maxgap
#     st      = np.zeros(ngap)
#     for j in range(ngap):
#         tx = []
#         for k in range(n2ave):
#             print('Using positions:')
#             print(k*jump,k*jump+gapvect[j])
#             tx.append(cu2[k*jump,:,:]-cu2[k*jump+gapvect[j],:,:])
#         st[j]=np.mean(np.mean(np.mean(tx,1),1))
#     return st



def comp_psd(img,  nbins=None, verbose=None):
    sx = (np.shape(img))[0]

    if nbins is None:
        nbins = sx//2
    img = img-np.mean(img)
    mask = np.invert(img.mask)
    img[mask ==0]=0
    tf2d = scipy.fft.fft2(img,norm='ortho')
    tf2d[0,0]=0
    tf2d_power_spectrum = np.abs(tf2d)**2

    kfreq = np.fft.fftfreq(sx) * sx                #freq spaziale in pixel
    kfreq2D = np.meshgrid(kfreq, kfreq)            #griglia di frequenze xx e yy
    knrm = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)  #griglia di frequenze distanza

    knrm = knrm.flatten()
    fourier_amplitudes = tf2d_power_spectrum.flatten()
    kbins = np.arange(0., nbins+1, 1.)
    kvals = 0.5 * (kbins[1:] + kbins[:-1])
    Abins, x_edges, y_edges = stats.binned_statistic(knrm, fourier_amplitudes,
                                         statistic = "sum",
                                         bins = nbins)

    #assert(Abins[0] == 0)
    #ediff = (np.sum(img[mask]i**2) - np.sum(Abins))/np.sum(img[mask]**2)
    ediff = (np.sum(img**2) - np.sum(Abins))/np.sum(img**2)

    if verbose is None:
        pass
    else:
        print("Energy difference %e" % ediff)
        print("RMS [nm] %5.2f" % (np.std(img)*1e9))
    Abins = Abins / np.sum(Abins) *(np.sum(img**2))
    return (kvals, Abins)

        


