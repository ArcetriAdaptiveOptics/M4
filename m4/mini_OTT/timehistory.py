import os
import glob
import numpy as np
import jdcal
from matplotlib.pyplot import *
import psutil
import scipy.fft
import scipy.stats as stats
from astropy.io import fits as pyfits
from m4.configuration import update_folder_paths as ufp
foldname = ufp.folders
from m4.ground import zernike, geo, read_data
from m4.ground.read_data import InterferometerConverter


ic = InterferometerConverter()
# a= foldname.BASE_PATH+'M4Data/OPTData/'  #'/mnt/data/M4/Data/M4Data/OPTData/'
a = foldname.OPT_DATA_FOLDER + "/"
# a='/mnt/m4storage/Data/M4Data/OPTData/'
# tn = '20210425_085242'   #fits
# tn  ='20210429_224400_noise'
from m4.configuration import read_4dconfig as readconf4d
from m4.ground.read_4DConfSettingFile import ConfSettingReader


# class TimeHist():
#     '''
#     Class to
#     HOW TO USE IT::
# f=
#
#     '''
#
# def __init__(self, tn):
#     """The constructor """
#     self.tracknum = tn
#
#     self._fold  = findTracknum(tn)
#     self._path = a+ self._fold
#     self._list = fileList(tn)
#     try:
#         self._confReader = ConfSettingReader(self.getConf4DSettingsPath)
#     finally:
#         pass
#
# def frame(self, id):
#
#     return frame(id, self._list)
#
# def averageFrames(self, start, stop):
#
#     return averageFrames(start, stop, self._list)
#
# def getConf4DSettingsPath(self):
#     fold = findTracknum(tn)
#     addfold = '/'
#     dirs = os.listdir(a+'/'+ fold+'/'+tn)
#     file_path = os.path.join(dirs, '4DSettings.ini')
#     return file_path


# @staticmethod
# def _storageFolder():
#    """ Creates the path where to save data"""
#    return fold_name.OPDSERIES


def findTracknum(tn):
    """
    Parameters
    ----------
    tn: string
        tracking number to be searched in the data folder

    Returns
    -------
    result: the specific data folder where the tracknum is found
    """

    # a= '/mnt/data/M4/Data/M4Data/OPTData/'
    lsdir = os.listdir(a)
    for i in lsdir:
        b = a + i
        z = os.listdir(b)
        check = False
        for j in z:
            check = j == tn
            if check == True:
                result = i
                return result

def _sortFunc4D(elem):
    iid = os.path.basename(elem)[:-3]
    iis = "%5.5i" % int(iid)
    return iis

def fileList(tn, fold=None, name=None):
    """
    Parameters
    ----------
    tn: str
        tracking number where to search for the images file list

    Returns
    -------
    lsdir : str | ArrayLike
        The list of image files
    """
    if fold is not None:
        if name is None:
            raise Exception
        tn = ""
        addfold = "/"
        # name = '*.4D'
        # addfold ='/hdf5/'
        fold1 = fold + "/" + tn + addfold
    else:

        fold = findTracknum(tn)
        addfold = "/"
        dirs = os.listdir(a + "/" + fold + "/" + tn)
        if dirs[0] == "hdf5":
            addfold = "/hdf5/"
            name = "img*.h5"
        elif dirs[0][-3:] == ".4D":
            name = "*.4D"
        else:
            name = "20*.fits"
        fold1 = a + "/" + fold + "/" + tn + addfold  # to be re-checked at OTT!!

    print(fold1 + name)
    lsdirs = glob.glob(fold1 + name)
    if name == "*.4D":
        lsdirs.sort(key=_sortFunc4D)
        lsdirs.sort(key=_sortFunc4D)
    else:
        lsdirs.sort()
    return lsdirs


def getConf4DSettingsPath(tn):
    fold = findTracknum(tn)
    addfold = "/"
    dirs = a + "/" + fold + "/" + tn
    file_path = dirs + "/4DSettings.ini"
    return file_path


def getCameraSettings(tn):
    file_path = getConf4DSettingsPath(tn)
    setting_reader = ConfSettingReader(file_path)
    width_pixel = setting_reader.getImageWidhtInPixels()
    height_pixel = setting_reader.getImageHeightInPixels()
    offset_x = setting_reader.getOffsetX()
    offset_y = setting_reader.getOffsetY()
    return [width_pixel, height_pixel, offset_x, offset_y]


def getFrameRate(tn):
    file_path = getConf4DSettingsPath(tn)
    setting_reader = ConfSettingReader(file_path)
    frame_rate = setting_reader.getFrameRate()
    return frame_rate


def read_phasemap(filename, thefold=None):

    # if thefold is not None:
    #    the = filename.split(thefold)[1]

    thetype = filename.split(".")[1]
    if thetype == "fits":
        hduList = pyfits.open(filename)
        img = hduList[0].data
        mask = hduList[1].data.astype(bool)

        """
        with pyfits.open(filename) as hduList:
        
            img = hduList[0].data
            mask = hduList[1].data
        """
        # mask = np.zeros(img.shape, dtype=bool)
        # mask[np.where(img == img.max())] = True

        img = np.ma.masked_array(img, mask)
        # hduList.close()

    if thetype == "4D":
        img = ic.fromPhaseCam6110(filename)

    if thetype == "h5":
        img = ic.fromPhaseCam4020(filename)

    return img


def averageFrames(first, last, fileList, thresh=None, fsel=None):
    """
    Parameters
    ----------
    first: first item
        tracking number where to search for the images file list

    Returns
    -------
    lsdir: the list of image files
    fold1: the path where the file is found
    """
    if fsel is None:
        fsel = np.arange(first, last + 1)

    # imcube = cubeFromList(fileList[x for x in flist])
    imcube = cubeFromList([fileList[x] for x in fsel])
    if thresh is None:
        aveimg = np.ma.mean(imcube, axis=0)

    else:
        img = imcube[0].data * 0
        mmask = imcube[0].mask
        mysize = imcube[0].compressed()
        nn = 0
        for i in imcube:
            if i.data.compressed.size > 1:
                nn += 1
                img += i.data
                mmask = np.ma.mask_or(i.mask, mmask)

        img = img / nn
        image = np.ma.masked_array(img, mask=mmask)

    return aveimg

def saveAverage(tn, id0=0, id1=None):
    fold = th.findTracknum(tn)
    fname = foldname.OPT_DATA_FOLDER + "/" + fold + "/" + tn + "/average.fits"
    print(fname)
    # check file esiste

    if os.path.isfile(fname) == True:
        print("average already exists")
    else:
        print("average do not exists jet")
        fl = fileList(tn)
        if id1 == None:
            id1 = np.size(fl) - 1
            # print(id0,id1)
        aveimg = averageFrames(id0, id1, fl)
        fold = findTracknum(tn)
        fname = foldname.OPT_DATA_FOLDER + "/" + fold + "/" + tn + "/average.fits"
        print(fname)
        pyfits.writeto(fname, aveimg.data)
        pyfits.append(fname, aveimg.mask.astype(np.uint8))


def openAverage(tn):
    fold = findTracknum(tn)
    fname = foldname.OPT_DATA_FOLDER + "/" + fold + "/" + tn + "/average.fits"

    if os.path.isfile(fname) == True:
        print("average exist")
        hduList = pyfits.open(fname)
        image = np.ma.masked_array(hduList[0].data, mask=hduList[1].data.astype(bool))

    else:
        print("average do not exist")

    return image


def removeZernike(ima, modes=np.array([1, 2, 3, 4])):

    coeff, mat = zernike.zernikeFit(ima, modes)
    surf = zernike.zernikeSurface(ima, coeff, mat)
    new_ima = ima - surf
    return new_ima


def removeZernikeAuxMask(img, mm, zlist):
    coef, mat = zernike.zernikeFitAuxmask(img, mm, zlist)
    surf = zernike.zernikeSurface(img, coef, mat)
    new_ima = np.ma.masked_array(img - surf, img.mask)
    return new_ima


def zernikePlot(mylist, modes=np.array(range(1, 11))):
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
    llist = fileList(tn)
    nfile = len(llist)
    npoints = int(nfile / gap) - 2
    slist = []
    for i in range(0, npoints):
        # print(i*gap)
        # print(i*gap+1)
        q0 = frame(i * gap, llist)
        q1 = frame(i * gap + 1, llist)
        diff = q1 - q0
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
    if np.size(nsig) == 1:
        thedim = 0

    else:
        thedim = 1

    if np.size(nsig) == 1:
        spe = np.fft.rfft(signal, norm="ortho")
        nn = np.sqrt(spe.shape[thedim])  # modRB

    else:
        spe = np.fft.rfft(signal, axis=1, norm="ortho")
        nn = np.sqrt(spe.shape[thedim])  # modRB
    spe = (np.abs(spe)) / nn
    freq = np.fft.rfftfreq(signal.shape[thedim], d=dt)

    if np.size(nsig) == 1:
        spe[0] = 0
    else:
        spe[:, 0] = 0
    if show is not None:
        for i in range(0, nn):
            plot(freq, spe[i, :])
    return spe, freq


def cubeFromList(fileList):
    image_list = []
    for i in fileList:
        ima = read_phasemap(i)
        image_list.append(ima)

    image_list = np.ma.masked_array(image_list)
    # print(psutil.Process().open_files())
    return image_list


def frame2ottFrame(img, croppar, flipOffset=True):
    off = croppar.copy()
    if flipOffset == True:
        off = np.flip(croppar)
        print("Offset values flipped:" + str(off))
    # nfullpix = np.array([2000,2000])
    nfullpix = np.array([2048, 2048])
    fullimg = np.zeros(nfullpix)
    fullmask = np.ones(nfullpix)
    offx = off[0]
    offy = off[1]
    sx = np.shape(img)[0]  # croppar[2]
    sy = np.shape(img)[1]  # croppar[3]
    fullimg[offx : offx + sx, offy : offy + sy] = img.data
    fullmask[offx : offx + sx, offy : offy + sy] = img.mask
    fullimg = np.ma.masked_array(fullimg, fullmask)
    return fullimg


def timevec(tn):
    fold = findTracknum(tn)
    flist = fileList(tn)
    nfile = len(flist)
    if fold == "OPDImages":
        tspace = 1.0 / 28.57
        timevec = range(nfile) * tspace

    if fold == "OPD_series":
        timevec = []
        for i in flist:
            pp = i.split(".")[0]
            tni = pp.split("/")[-1]
            y = tni[0:4]
            mo = tni[4:6]
            d = tni[6:8]
            h = float(tni[9:11])
            mi = float(tni[11:13])
            s = float(tni[13:15])
            jdi = sum(jdcal.gcal2jd(y, mo, d)) + h / 24 + mi / 1440 + s / 86400
            timevec.append(jdi)
        timevec = np.array(timevec)

    return timevec


def track2jd(tni):
    y, mo, d, h, mi, s = track2date(tni)
    jdi = sum(jdcal.gcal2jd(y, mo, d)) + h / 24 + mi / 1440 + s / 86400
    return jdi


def track2date(tni):
    y = tni[0:4]
    mo = tni[4:6]
    d = tni[6:8]
    h = float(tni[9:11])
    mi = float(tni[11:13])
    s = float(tni[13:15])
    return y, mo, d, h, mi, s


def runningMean(vec, npoints):

    return np.convolve(vec, np.ones(npoints), "valid") / npoints


def readTemperatures(tn):
    fold = findTracknum(tn)
    fname = foldname.OPT_DATA_FOLDER + "/" + fold + "/" + tn + "/temperature.fits"
    temperatures = (pyfits.open(fname))[0].data
    return temperatures


def readZernike(tn):
    fold = findTracknum(tn)
    fname = foldname.OPT_DATA_FOLDER + "/" + fold + "/" + tn + "/zernike.fits"
    temperatures = (pyfits.open(fname))[0].data
    return temperatures


def strfunct(vect, gapvect):
    """
    vect shall be npoints x m
    the strfunct is calculate m times over the npoints time series
    returns stf(n_timeseries x ngaps)
    """
    nn = np.shape(vect)
    maxgap = np.max(gapvect)
    ngap = len(gapvect)
    n2ave = int(nn / (maxgap)) - 1  # or -maxgap??
    jump = maxgap
    st = np.zeros(ngap)
    for j in range(ngap):
        tx = []
        for k in range(n2ave):
            print("Using positions:")
            print(k * jump, k * jump + gapvect[j])
            tx.append((vect[k * jump] - vect[k * jump + gapvect[j]]) ** 2)
        st[j] = np.mean(np.sqrt(tx))
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


def readFrameCrop(tn):
    fold = findTracknum(tn)
    path = a + "/" + fold + "/" + tn
    print(path)
    w, h, x, y, fr = readconf4d.getCameraConfig(path)
    return np.array([w, h, x, y], dtype=np.int64)  #modRB20240518, was dtype=int64



def readFrameRate(tn):
    fold = findTracknum(tn)
    path = a + "/" + fold + "/" + tn
    print(path)
    w, h, x, y, fr = readconf4d.getCameraConfig(tn)
    return fr


def comp_filtered_image(
    imgin, verbose=False, disp=False, d=1, crop=True, freq2filter=None
):

    #    if crop:
    #        cir = geo.qpupil(-1*imgin.mask+1)
    #        cir = np.array(cir[0:3]).astype(int)
    #        img = imgin.data[cir[0]-cir[2]:cir[0]+cir[2],cir[1]-cir[2]:cir[1]+cir[2]]
    #        m = imgin.mask[cir[0]-cir[2]:cir[0]+cir[2],cir[1]-cir[2]:cir[1]+cir[2]]
    #        img = np.ma.masked_array(img, m)
    #    else:
    img = imgin.copy()

    sx = (np.shape(img))[0]

    # img = img-np.mean(img)
    mask = np.invert(img.mask)

    img[mask == 0] = 0
    norm = "ortho"
    tf2d = scipy.fft.fft2(img.data, norm=norm)
    # tf2d[0,0]=0

    kfreq = np.fft.fftfreq(sx, d=d)  # frequenza in cicli
    kfreq2D = np.meshgrid(kfreq, kfreq)  # griglia di frequenze xx e yy
    knrm = np.sqrt(kfreq2D[0] ** 2 + kfreq2D[1] ** 2)  # griglia di frequenze distanza

    # maschera (da mettere opzionale) per limitarsi al cerchio e non prendere il quadrato
    fmask1 = 1.0 * (knrm > np.max(kfreq))

    if freq2filter is None:
        fmin = -1
        fmax = np.max(kfreq)
    else:
        fmin, fmax = freq2filter

    fmask2 = 1.0 * (knrm > fmax)
    fmask3 = 1.0 * (knrm < fmin)

    fmask = (fmask1 + fmask2 + fmask3) > 0

    tf2d_filtered = tf2d.copy()
    tf2d_filtered[fmask] = 0
    imgf = scipy.fft.ifft2(tf2d_filtered, norm=norm)
    imgout = np.ma.masked_array(np.real(imgf), mask=imgin.mask)
    if disp:
        figure()
        imshow(knrm)
        title("freq")
        figure()
        imshow(fmask1)
        title("fmask1")
        figure()
        imshow(fmask2)
        title("fmask2")
        figure()
        imshow(fmask3)
        title("fmask3")
        figure()
        imshow(fmask)
        title("fmask")
        figure()
        imshow(np.abs(tf2d))
        title("Initial spectrum")
        figure()
        imshow(np.abs(tf2d_filtered))
        title("Filtered spectrum")
        figure()
        imshow(imgin)
        title("Initial image")
        figure()
        imshow(imgout)
        title("Filtered image")

    e1 = np.sqrt(np.sum(img[mask] ** 2) / np.sum(mask)) * 1e9
    e2 = np.sqrt(np.sum(imgout[mask] ** 2) / np.sum(mask)) * 1e9
    e3 = np.sqrt(np.sum(np.abs(tf2d) ** 2) / np.sum(mask)) * 1e9
    e4 = np.sqrt(np.sum(np.abs(tf2d_filtered) ** 2) / np.sum(mask)) * 1e9

    if verbose:
        print("RMS image [nm]            %5.2f" % e1)
        print("RMS image filtered [nm]   %5.2f" % e2)
        print("RMS spectrum              %5.2f" % e3)
        print("RMS spectrum filtered     %5.2f" % e4)

    return imgout


def comp_psd(
    imgin,
    nbins=None,
    norm="backward",
    verbose=False,
    disp=False,
    d=1,
    sigma=None,
    crop=True,
):

    if crop:
        cir = geo.qpupil(-1 * imgin.mask + 1)
        cir = np.array(cir[0:3]).astype(int)
        img = imgin.data[
            cir[0] - cir[2] : cir[0] + cir[2], cir[1] - cir[2] : cir[1] + cir[2]
        ]
        m = imgin.mask[
            cir[0] - cir[2] : cir[0] + cir[2], cir[1] - cir[2] : cir[1] + cir[2]
        ]
        img = np.ma.masked_array(img, m)
    else:
        img = imgin.copy()

    sx = (np.shape(img))[0]

    if nbins is None:
        nbins = sx // 2
    img = img - np.mean(img)
    mask = np.invert(img.mask)

    # DA IMPLEMENTARE: finestre contro effetti di bordo
    #        maskf = scipy.ndimage.fourier_gaussian(np.array(mask, dtype=float), sigma=sigma)

    #    img=maskf*img#0
    img[mask == 0] = 0
    if sigma is not None:
        img = scipy.ndimage.fourier_gaussian(img, sigma=sigma)
    tf2d = scipy.fft.fft2(img, norm=norm)
    tf2d[0, 0] = 0
    tf2d_power_spectrum = np.abs(tf2d) ** 2

    # kfreq = np.fft.fftfreq(sx) * sx                #freq spaziale in pixel
    kfreq = np.fft.fftfreq(sx, d=d)  # frequenza in cicli
    kfreq2D = np.meshgrid(kfreq, kfreq)  # griglia di frequenze xx e yy
    knrm = np.sqrt(kfreq2D[0] ** 2 + kfreq2D[1] ** 2)  # griglia di frequenze distanza

    # maschera (da mettere opzionale) per limitarsi al cerchio e non prendere il quadrato
    fmask = knrm < np.max(kfreq)
    # fmask = mask.copy()
    knrm = knrm[fmask].flatten()
    fourier_amplitudes = tf2d_power_spectrum[fmask].flatten()
    Abins, x_edges, y_edges = stats.binned_statistic(
        knrm, fourier_amplitudes, statistic="sum", bins=nbins
    )

    # assert(Abins[0] == 0)
    # ediff = (np.sum(img[mask]i**2) - np.sum(Abins))/np.sum(img[mask]**2)
    e1 = np.sum(img[mask] ** 2 / np.sum(mask))
    e2 = np.sum(Abins) / np.sum(mask)
    ediff = np.abs(e2 - e1) / e1

    fout = kfreq[0 : sx // 2]
    Aout = Abins / np.sum(mask)
    if verbose:
        print("RMS from spectrum %e" % np.sqrt(e2))
        print("RMS [nm]          %5.2f" % (np.std(img[mask]) * 1e9))
    else:
        print("Sampling          %e" % d)
        print("Energy signal     %e" % e1)
        print("Energy spectrum    %e" % e2)
        print("Energy difference %e" % ediff)
        print("RMS from spectrum %e" % np.sqrt(e2))
        print("RMS [nm]          %5.2f" % (np.std(img[mask]) * 1e9))
        print(kfreq[0:4])
        print(kfreq[-4:])
    if disp == True:
        figure()
        plot(fout[1:], Aout[1:] * fout[1:], ".")
        yscale("log")
        xscale("log")
        title("Power spectrum")
        xlabel("[Hz]")
        ylabel("[A^2]")

    return fout, Aout


def comp_psd_old(img, nbins=None, verbose=None):
    sx = (np.shape(img))[0]

    if nbins is None:
        nbins = sx // 2
    img = img - np.mean(img)
    mask = np.invert(img.mask)
    img[mask == 0] = 0
    tf2d = scipy.fft.fft2(img, norm="ortho")
    tf2d[0, 0] = 0
    tf2d_power_spectrum = np.abs(tf2d) ** 2

    kfreq = np.fft.fftfreq(sx) * sx  # freq spaziale in pixel
    kfreq2D = np.meshgrid(kfreq, kfreq)  # griglia di frequenze xx e yy
    knrm = np.sqrt(kfreq2D[0] ** 2 + kfreq2D[1] ** 2)  # griglia di frequenze distanza

    knrm = knrm.flatten()
    fourier_amplitudes = tf2d_power_spectrum.flatten()
    kbins = np.arange(0.0, nbins + 1, 1.0)
    kvals = 0.5 * (kbins[1:] + kbins[:-1])
    Abins, x_edges, y_edges = stats.binned_statistic(
        knrm, fourier_amplitudes, statistic="sum", bins=nbins
    )

    # assert(Abins[0] == 0)
    # ediff = (np.sum(img[mask]i**2) - np.sum(Abins))/np.sum(img[mask]**2)
    ediff = (np.sum(img**2) - np.sum(Abins)) / np.sum(img**2)

    if verbose is None:
        pass
    else:
        print("Energy difference %e" % ediff)
        print("RMS [nm] %5.2f" % (np.std(img) * 1e9))
    Abins = Abins / np.sum(Abins) * (np.sum(img**2))
    return (kvals, Abins)


def integrate_psd(y, img):
    nn = np.sqrt(np.sum(-1 * img.mask + 1))
    yint = np.sqrt(np.cumsum(y)) / nn
    return yint


def tnscan(tn0, tn1):
    """
        returns the list of tracknum in a given datafolder, in between tn0 and tn1
    syntax: tnlist = tnscan(tn0, tn1)
    """
    datafold = findTracknum(tn0)
    ll = sorted(os.listdir(foldname.OPT_DATA_FOLDER + "/" + datafold))
    id0 = ll.index(tn0)
    id1 = ll.index(tn1)
    tnlist = ll[id0 : id1 + 1]
    return tnlist
