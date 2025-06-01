"""
Author(s)
---------
    - Runa Briguglio: created 2020
    - Pietro Ferraiuolo: modified 2024

Description
-----------

"""

import os
import jdcal
import numpy as np
from scipy import stats, fft, ndimage
import matplotlib.pyplot as plt
from m4.utils import osutils as osu
from m4.configuration import update_folder_paths as ufp
from m4.ground import zernike, geo, read_data as rd
from m4.ground.read_data import InterferometerConverter

foldname = ufp.folders
ic = InterferometerConverter()
OPTDATA = foldname.OPT_DATA_FOLDER
OPDIMG = foldname.OPD_IMAGES_ROOT_FOLDER
OPDSER = foldname.OPD_SERIES_ROOT_FOLDER


def averageFrames(
    tn: str,
    first: int = None,
    last: int = None,
    file_selector: list = None,
    thresh: bool = False,
):
    """
    Perform the average of a list of images, retrievable through a tracking
    number.

    Parameters
    ----------
    tn : str
        Data Tracking Number.
    first : int, optional
        Index number of the first file to consider. If None, the first file in
        the list is considered.
    last : int, optional
        Index number of the last file to consider. If None, the last file in
        list is considered.
    file_selector : list, optional
        A list of integers, representing the specific files to load. If None,
        the range (first->last) is considered.
    thresh : bool, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    aveimg : ndarray
        Final image of averaged frames.

    """
    fileList = osu.getFileList(tn, fold=OPDSER, key="20")
    if first is not None and last is not None:
        fl = [
            fileList[x]
            for x in np.arange(first, last, 1)
            if file_selector is None or x in file_selector
        ]
    else:
        first = 0
        last = len(fileList)
        fl = [
            fileList[x]
            for x in np.arange(first, last, 1)
            if file_selector is None or x in file_selector
        ]
    imcube = osu.createCube(fl)
    if thresh is False:
        aveimg = np.ma.mean(imcube, axis=2)
    else:
        img = imcube[:, :, 0].data * 0
        mmask = imcube[:, :, 0].mask
        nn = 0
        for j in range(imcube.shape[2]):
            im = imcube[:, :, j]
            size = im.data.compressed.size
            if size > 1:
                nn += 1
                img += im.data
                mmask = np.ma.mask_or(im.mask, mmask)
        img = img / nn
        aveimg = np.ma.masked_array(img, mask=mmask)
    return aveimg


def saveAverage(tn, average_img=None, overwrite: bool = False, **kwargs):
    """
    Saves an averaged frame, in the same folder as the original frames. If no
    averaged image is passed as argument, it will create a new average for the
    specified tracking number, and additional arguments, the same as ''averageFrames''
    can be specified.

    Parameters
    ----------
    tn : str
        Tracking number where to save the average frame file. If average_img is
        None, it is the tracking number of the data that will be averaged
    average_img : ndarray, optional
        Result average image of multiple frames. If it's None, it will be generated
        from data found in the tracking number folder. Additional arguments can
        be passed on
    **kwargs : additional optional arguments
        The same arguments as ''averageFrames'', to specify the averaging method.

        tn : str
            Data Tracking Number.
        first : int, optional
            Index number of the first file to consider. If None, the first file in
            the list is considered.
        last : int, optional
            Index number of the last file to consider. If None, the last file in
            list is considered.
        file_selector : list, optional
            A list of integers, representing the specific files to load. If None,
            the range (first->last) is considered.
        thresh : bool, optional
            DESCRIPTION. The default is None.
    """
    fname = os.path.join(OPDSER, tn, "average.fits")
    if os.path.isfile(fname):
        print(f"Average '{fname}' already exists")
    else:
        if average_img is None:
            first = kwargs.get("first", None)
            last = kwargs.get("last", None)
            fsel = kwargs.get("file_selector", None)
            thresh = kwargs.get("tresh", False)
            average_img = averageFrames(
                tn, first=first, last=last, file_selector=fsel, thresh=thresh
            )
        rd.save_phasemap(fname, average_img, overwrite=overwrite)
        print(f"Saved average at '{fname}'")


def openAverage(tn):
    """
    Loads an averaged frame from an 'average.fits' file, found inside the input
    tracking number

    Parameters
    ----------
    tn : str
        Tracking number of the averaged frame.

    Returns
    -------
    image : ndarray
        Averaged image.

    Raises
    ------
    FileNotFoundError
        Raised if the file does not exist.
    """
    fname = os.path.join(OPDSER, tn, "average.fits")
    try:
        image = rd.readFits_maskedImage(fname)
        print(f"Average loaded: '{fname}'")
    except FileNotFoundError as err:
        raise FileNotFoundError(f"Average file '{fname}' does not exist!") from err
    return image


def runningDiff(tn, gap=2):
    """


    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.
    gap : TYPE, optional
        DESCRIPTION. The default is 2.

    Returns
    -------
    svec : TYPE
        DESCRIPTION.

    """
    llist = osu.getFileList(tn)
    nfile = len(llist)
    npoints = int(nfile / gap) - 2
    slist = []
    for i in range(0, npoints):
        q0 = frame(i * gap, llist)
        q1 = frame(i * gap + 1, llist)
        diff = q1 - q0
        diff = zernike.removeZernike(diff)
        slist.append(diff.std())
    svec = np.array(slist)
    return svec


def frame(idx, mylist):
    """


    Parameters
    ----------
    id : TYPE
        DESCRIPTION.
    mylist : TYPE
        DESCRIPTION.

    Returns
    -------
    img : TYPE
        DESCRIPTION.

    """
    mytype = type(mylist)
    if mytype is list:
        img = rd.read_phasemap(mylist[idx])
    if mytype is np.ma.core.MaskedArray:
        img = mylist[idx]
    return img


def spectrum(signal, dt=1, show=None):
    """


    Parameters
    ----------
    signal : ndarray
        DESCRIPTION.
    dt : float, optional
        DESCRIPTION. The default is 1.
    show : bool, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    spe : float | ndarray
        DESCRIPTION.
    freq : float | ArrayLike
        DESCRIPTION.

    """
    # https://numpy.org/doc/stable/reference/generated/numpy.angle - Spectrum phase
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
        plt.figure()
        for i in range(0, len(spe)):
            plt.plot(freq, spe[i, :], label=f"Channel {i}")
        plt.xlabel(r"Frequency [$Hz$]")
        plt.ylabel("PS Amplitude")
        plt.legend(loc="best")
        plt.show()
    return spe, freq


def frame2ottFrame(img, croppar, flipOffset=True):
    """


    Parameters
    ----------
    img : TYPE
        DESCRIPTION.
    croppar : TYPE
        DESCRIPTION.
    flipOffset : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    fullimg : TYPE
        DESCRIPTION.

    """
    off = croppar.copy()
    if flipOffset is True:
        off = np.flip(croppar)
        print(f"Offset values flipped: {str(off)}")
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
    """


    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.

    Returns
    -------
    timevector : TYPE
        DESCRIPTION.

    """
    fold = osu.findTracknum(tn)
    flist = osu.getFileList(tn)
    nfile = len(flist)
    if "OPDImages" in fold:
        tspace = 1.0 / 28.57
        timevector = range(nfile) * tspace
    elif "OPD_series" in fold:
        timevector = []
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
            timevector.append(jdi)
        timevector = np.array(timevec)
    return timevector


def track2jd(tni):
    """


    Parameters
    ----------
    tni : TYPE
        DESCRIPTION.

    Returns
    -------
    jdi : TYPE
        DESCRIPTION.

    """
    t = track2date(tni)
    jdi = sum(jdcal.gcal2jd(t[0], t[1], t[2])) + t[3] / 24 + t[4] / 1440 + t[5] / 86400
    return jdi


def track2date(tni):
    """
    Converts a tracing number into a list containing year, month, day, hour,
    minutes and seconds, divied.

    Parameters
    ----------
    tni : str
        Tracking number to be converted.

    Returns
    -------
    time : list
        List containing the date element by element.
        [0] y : str
            Year.
        [1] mo : str
            Month.
        [2] d : str
            Day.
        [3] h : float
            Hour.
        [4] mi : float
            Minutes.
        [5] s : float
            Seconds.
    """
    y = tni[0:4]
    mo = tni[4:6]
    d = tni[6:8]
    h = float(tni[9:11])
    mi = float(tni[11:13])
    s = float(tni[13:15])
    time = [y, mo, d, h, mi, s]
    return time


def runningMean(vec, npoints):
    """


    Parameters
    ----------
    vec : TYPE
        DESCRIPTION.
    npoints : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return np.convolve(vec, np.ones(npoints), "valid") / npoints


def readTemperatures(tn):
    """


    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.

    Returns
    -------
    temperatures : TYPE
        DESCRIPTION.

    """
    fold = osu.findTracknum(tn, complete_path=True)
    fname = os.path.join(fold, tn, "temperature.fits")
    temperatures = rd.readFits_data(fname)
    return temperatures


def readZernike(tn):
    """


    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.

    Returns
    -------
    temperatures : TYPE
        DESCRIPTION.

    """
    fold = osu.findTracknum(tn, complete_path=True)
    fname = os.path.join(fold, tn, "zernike.fits")
    zernikes = rd.readFits_data(fname)
    return zernikes


def zernikePlot(mylist, modes=np.array(range(1, 11))):
    """


    Parameters
    ----------
    mylist : TYPE
        DESCRIPTION.
    modes : TYPE, optional
        DESCRIPTION. The default is np.array(range(1, 11)).

    Returns
    -------
    zcoeff : TYPE
        DESCRIPTION.

    """
    mytype = type(mylist)
    if mytype is list:
        imgcube = osu.createCube(mylist)
    if mytype is np.ma.core.MaskedArray:
        imgcube = mylist
    zlist = []
    for i in range(len(imgcube)):
        print(i)
        coeff, _ = zernike.zernikeFit(imgcube[i], modes)
        zlist.append(coeff)
    zcoeff = np.array(zlist)
    zcoeff = zcoeff.T
    return zcoeff


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


def comp_filtered_image(imgin, verbose=False, disp=False, d=1, freq2filter=None):
    """


    Parameters
    ----------
    imgin : TYPE
        DESCRIPTION.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.
    disp : TYPE, optional
        DESCRIPTION. The default is False.
    d : TYPE, optional
        DESCRIPTION. The default is 1.
    freq2filter : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    imgout : TYPE
        DESCRIPTION.
    """
    img = imgin.copy()
    sx = (np.shape(img))[0]
    mask = np.invert(img.mask)
    img[mask == 0] = 0
    norm = "ortho"
    tf2d = fft.fft2(img.data, norm=norm)
    kfreq = np.fft.fftfreq(sx, d=d)  # frequency in cicles
    kfreq2D = np.meshgrid(kfreq, kfreq)  # frequency grid x,y
    knrm = np.sqrt(kfreq2D[0] ** 2 + kfreq2D[1] ** 2)  # freq. grid distance
    # TODO optional mask to get the circle and not the square
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
    imgf = fft.ifft2(tf2d_filtered, norm=norm)
    imgout = np.ma.masked_array(np.real(imgf), mask=imgin.mask)
    if disp:
        plt.figure()
        plt.imshow(knrm)
        plt.title("freq")
        plt.figure()
        plt.imshow(fmask1)
        plt.title("fmask1")
        plt.figure()
        plt.imshow(fmask2)
        plt.title("fmask2")
        plt.figure()
        plt.imshow(fmask3)
        plt.title("fmask3")
        plt.figure()
        plt.imshow(fmask)
        plt.title("fmask")
        plt.figure()
        plt.imshow(np.abs(tf2d))
        plt.title("Initial spectrum")
        plt.figure()
        plt.imshow(np.abs(tf2d_filtered))
        plt.title("Filtered spectrum")
        plt.figure()
        plt.imshow(imgin)
        plt.title("Initial image")
        plt.figure()
        plt.imshow(imgout)
        plt.title("Filtered image")
    e1 = np.sqrt(np.sum(img[mask] ** 2) / np.sum(mask)) * 1e9
    e2 = np.sqrt(np.sum(imgout[mask] ** 2) / np.sum(mask)) * 1e9
    e3 = np.sqrt(np.sum(np.abs(tf2d) ** 2) / np.sum(mask)) * 1e9
    e4 = np.sqrt(np.sum(np.abs(tf2d_filtered) ** 2) / np.sum(mask)) * 1e9
    if verbose:
        print(f"RMS image [nm]            {e1:.2f}")
        print(f"RMS image filtered [nm]   {e2:.2f}")
        print(f"RMS spectrum              {e3:.2f}")
        print(f"RMS spectrum filtered     {e4:.2f}")
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
    """


    Parameters
    ----------
    imgin : TYPE
        DESCRIPTION.
    nbins : TYPE, optional
        DESCRIPTION. The default is None.
    norm : TYPE, optional
        DESCRIPTION. The default is "backward".
    verbose : TYPE, optional
        DESCRIPTION. The default is False.
    disp : TYPE, optional
        DESCRIPTION. The default is False.
    d : TYPE, optional
        DESCRIPTION. The default is 1.
    sigma : TYPE, optional
        DESCRIPTION. The default is None.
    crop : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    fout : TYPE
        DESCRIPTION.
    Aout : TYPE
        DESCRIPTION.

    """
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
    # TODO: finestre contro effetti di bordo
    #        maskf = scipy.ndimage.fourier_gaussian(np.array(mask, dtype=float), sigma=sigma)
    #    img=maskf*img#0
    img[mask == 0] = 0
    if sigma is not None:
        img = ndimage.fourier_gaussian(img, sigma=sigma)
    tf2d = fft.fft2(img, norm=norm)
    tf2d[0, 0] = 0
    tf2d_power_spectrum = np.abs(tf2d) ** 2
    kfreq = np.fft.fftfreq(sx, d=d)  # frequency in cicles
    kfreq2D = np.meshgrid(kfreq, kfreq)  # freq. grid
    knrm = np.sqrt(kfreq2D[0] ** 2 + kfreq2D[1] ** 2)  # freq. grid distance
    # TODO optional mask to get the circle and not the square
    fmask = knrm < np.max(kfreq)
    knrm = knrm[fmask].flatten()
    fourier_amplitudes = tf2d_power_spectrum[fmask].flatten()
    Abins, _, _ = stats.binned_statistic(
        knrm, fourier_amplitudes, statistic="sum", bins=nbins
    )
    e1 = np.sum(img[mask] ** 2 / np.sum(mask))
    e2 = np.sum(Abins) / np.sum(mask)
    ediff = np.abs(e2 - e1) / e1
    fout = kfreq[0 : sx // 2]
    Aout = Abins / np.sum(mask)
    if verbose:
        print(f"Sampling          {d:}")
        print(f"Energy signal     {e1}")
        print(f"Energy spectrum   {e2}")
        print(f"Energy difference {ediff}")
        print(f"RMS from spectrum {np.sqrt(e2)}")
        print(f"RMS [nm]          {(np.std(img[mask])*1e9):.2f}")
        print(kfreq[0:4])
        print(kfreq[-4:])
    else:
        print(f"RMS from spectrum {np.sqrt(e2)}")
        print(f"RMS [nm]          {(np.std(img[mask])*1e9):.2f}")
    if disp is True:
        plt.figure()
        plt.plot(fout[1:], Aout[1:] * fout[1:], ".")
        plt.yscale("log")
        plt.xscale("log")
        plt.title("Power spectrum")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("Amplitude [A^2]")
    return fout, Aout


def integrate_psd(y, img):
    nn = np.sqrt(np.sum(-1 * img.mask + 1))
    yint = np.sqrt(np.cumsum(y)) / nn
    return yint
