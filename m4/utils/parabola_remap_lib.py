import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from arte.utils import rebin
import opticalib
import opticalib.analyzer as th # from m4.analyzers import timehistory as th
# from m4.ground import geo, zernike as zern, read_data as rd
from opticalib.ground import geometry as geo # from m4.ground import geo
from opticalib.ground import modal_decomposer as md
# from opticalib.ground.osutils import newtn
from opticalib.ground import osutils as osu # from m4.utils import osutils as osu

#from m4.configuration import folders as foldname  modRB 20261008 replaced with opt.folders
foldname = opticalib.folders 
from m4.ground import read_ottcalib_conf as roc
import m4.utils.parabola_footprint_registration as pr
from m4.utils.parabola_identification import ParabolaActivities
from m4 import userscripts as usr  #temporary, we need intoFullFrame
zern = md.ZernikeFitter() # this is a guess, MM260121

pa = ParabolaActivities()
pfr = pr.ParabolaFootprintRegistration()
OPDSERIES = foldname.OPD_SERIES_ROOT_FOLDER

zern2remove_after = [1,2,3,4]
marksizecgh = 24
marksizeott = 28


def register_parabola2ott(tnconf, show=False, forder=10, dosave = False):
    cgh_image, ott_image, cghf, ottf,filtinfo = init_data(tnconf)
    if show is not False:
        view_markers(cghf, ottf)
    cgh_tra = par_remap(cgh_image,  cghf, ottf, forder=forder)
    cgh_tra = np.ma.masked_array(cgh_tra.data, cgh_tra == 0)
    cgh_tra = zern.remove_zernike(cgh_tra, zern2remove_after)
    #par_filtered = az.comp_filtered_image(par_remapped,  d=filtinfo[2], verbose=True, disp=False, freq2filter=filtinfo[0:2])

    tn = 'not saved'
    if dosave is True:
        tn = save_registration(cgh_tra, tnconf)
    print('Saved a new registration in '+tn)
    return cgh_tra

def par_remap(cgh_image, cghf, ottf, forder=10):
    _, _, cgh_tra, _ = pfr.image_transformation(
        cgh_image, cgh_image * 0, cghf, ottf, forder=forder
    )
    cgh_tra = np.ma.masked_array(cgh_tra.data, cgh_tra == 0)
    return cgh_tra

def save_registration(img, tnconf):  # (img,cgh_tn_img,cgh_tn_marker,ott_tn_marker):
    tn = osu.newtn()
    print(tn)
    fold = foldname.PARABOLA_REMAPPED_FOLDER + "/" + tn + "/"
    os.mkdir(fold)
    name = fold + "par_remapped.fits"
    pyfits.writeto(name, img.data)
    pyfits.append(name, img.mask.astype(int))
    copyConf(tnconf, tn)
    return tn

def init_data(tnconf):
    (
        cgh_tn_marker,
        cgh_tn_img,
        tnpar,
        mark_cgh_list,
        f0,
        f1,
        ott_tn_marker,
        ott_tn_img,
        mark_ott_list,
        px_ott,
    ) = read_init_data(tnconf)
    cghf = marker_data(cgh_tn_marker, mark_cgh_list, marksizecgh, flip=True)
    # here modified
    ntn = len(ott_tn_marker)
    if ntn > 1:
        print("Multi Tracknum")
        nlen = []
        for i in mark_ott_list:
            nlen.append(len(i))
        print("N tracknum" + str(ntn))
        ottf = np.zeros([2, np.sum(nlen)])
        pos = 0
        for i in np.arange(ntn):
            ottf[:, pos : pos + nlen[i]] = marker_data(
                ott_tn_marker[i], mark_ott_list[i], marksizeott, flip=False
            )
            pos = pos + nlen[i]
    else:
        ottf = marker_data(ott_tn_marker[0], mark_ott_list, marksizeott, flip=False)
    cgh_image = image_data(cgh_tn_img, flip=True)
    ott_image = image_data(ott_tn_img, flip=False)
    filtering_info = [f0, f1, px_ott]
    return cgh_image, ott_image, cghf, ottf, filtering_info


def marker_data(tn_marker, mark_list, diam, flip=False):
    """
    usage:
    for ott markers: marker_data(tn,mark_list, 28,flip=False)
    for cgh markers: marker_data(tn,mark_list, 24,flip=True)
    if tn_marker is a tnvector, mark list shall be a 2D vector
    """
    """for ott markers: marker_data(tn,mark_list, 28,flip=False)
    for cgh markers: marker_data(tn,mark_list, 24,flip=True)
    if tn_marker is a tnvector, mark list shall be a 2D vector')"""
    fl0 = osu.get_file_list(tn_marker, key='20')
    img0 = opticalib.read_phasemap(fl0[0])  #th.frame(0, fl0)
    if flip is True:
        img0 = np.fliplr(img0)
        print("flipping the frame")
    off_marker = (opticalib.get_camera_settings(tn_marker))[2:4]

    p0 = getMarkers(tn_marker, flip, diam)
    p0 = coord2ottcoord(p0, off_marker)
    if mark_list is not None:
        p0 = p0[:, mark_list]
    return p0

def image_data(tn_img, flip=False):
    path = foldname.OPD_SERIES_ROOT_FOLDER + "/"
    fname = path + tn_img + "/average.fits"
    img = opticalib.read_phasemap(fname)
    if flip is True:
        img = np.fliplr(img)
    conf = opticalib.get_camera_settings(tn_img)
    offs = conf[2:4]
    img = usr.into_full_frame(img, offs)
    return img


def view_markers(p0, p1):
    fig, ax = plt.subplots()  # figure()
    plt.plot(p0[1, :], p0[0, :], "o")
    ax.axis("equal")
    plt.plot(p1[1, :], p1[0, :], "x")
    for i in range(np.shape(p0)[1]):
        ax.text(
            p0[1, i], p0[0, i], str(i), color='b'  ) 
    for i in range(np.shape(p1)[1]):  
        ax.text(
            p1[1, i] + 10, p1[0, i] + 10, str(i), color='r'  ) 
    plt.title("Markers position comparison")
    plt.show()
    p00 = marker_remap(p0, p1)
    dd = np.sqrt((p00[0, :] - p1[0, :]) ** 2 + (p00[1, :] - p1[1, :]) ** 2)
    fig, ax, plt.subplots()  # figure()
    plt.scatter(p1[0, :], p1[1, :], dd * 50, dd)
    for i in range(np.shape(p0)[1]):
        ax.text(p1[0, i], p0[1, i], str(i))
    plt.title("Remapping error")
    plt.colorbar()


def marker_remap(cghf, ottf, forder=10):
    polycoeff = pfr.fit_trasformation_parameter(cghf, ottf,forder)
    base_cgh = pfr._expandbase(cghf[0, :], cghf[1, :])
    cghf_tra = np.transpose(np.dot(np.transpose(base_cgh), np.transpose(polycoeff)))
    return cghf_tra



def getMarkers(tn, flip=False, diam=24, thr=0.2):
    npix = 3.14 * (diam / 2) ** 2
    fl = osu.get_file_list(tn, key='20')
    nf = len(fl)
    pos = np.zeros([2, 25, nf])
    for j in range(nf):
        img = th.frame(j, fl)
        if flip is True:
            print("flipping")
            img = np.fliplr(img)
        imaf = pa.rawMarkersPos(img)
        c0 = pa.filterMarkersPos(imaf, (1 - thr) * npix, (1 + thr) * npix)
        if j == 0:
            nmark = np.shape(c0)[1]
            pos = np.zeros([2, nmark, nf])
        pos[:, :, j] = c0
    pos = np.average(pos, 2)
    return pos


def coord2ottcoord(vec1, off, flipOffset=True):
    off1 = off.copy()
    if flipOffset == True:
        off1 = np.flip(off)
        print("Offset values flipped:" + str(off1))
    vec = vec1.copy()
    for ii in range(np.shape(vec)[1]):
        vec[:, ii] = vec[:, ii] + off1
    return vec




def marker_general_remap(cghf, ottf, pos2t):
    """
    transforms the pos2t coordinates, using the cghf and ottf coordinates to create the trnasformation
    """
    polycoeff = pfr.fit_trasformation_parameter(cghf, ottf)
    base_cgh = pfr._expandbase(pos2t[0, :], pos2t[1, :])
    cghf_tra = np.transpose(np.dot(np.transpose(base_cgh), np.transpose(polycoeff)))
    return cghf_tra


def marker_data_all(
    cgh_tn_marker, ott_tn_marker, off_cgh_marker, off_ott_marker, mark_cgh, mark_ott
):
    tn0 = cgh_tn_marker
    tn1 = ott_tn_marker
    fl0 = osu.get_file_list(tn0, key='20')
    img0 = opticalib.read_phasemap(fl0[0]) 
    img0 = np.fliplr(img0)

    p0 = getMarkers(tn0, flip=True, diam=marksizecgh)
    p1 = getMarkers(tn1, diam=marksizeott)

    p0 = coord2ottcoord(p0, off_cgh_marker)
    p1 = coord2ottcoord(p1, off_ott_marker)
    pcgh = p0[:, mark_cgh]
    pott = p1[:, mark_ott]
    return pcgh, pott



def plot_markers(p0):
    _, ax = plt.subplots()
    plt.plot(p0[1, :], p0[0, :], "o")
    ax.axis("equal")
    plt.title("Markers position in the frame")
    for i in range(np.shape(p0)[1]):
        ax.text(p0[1, i], p0[0, i], str(i))


def markers_explorer(tn):
    p0 = marker_data(tn, None, diam=28, flip=False)
    plot_markers(p0)
    plt.xlim(0, 2048)
    plt.ylim(0, 2048)
    plt.title(tn)



def read_init_data(tnconf):
    (
        cgh_tn_marker,
        cgh_tn_img,
        tnpar,
        mark_cgh_list,
        f0,
        f1,
        ott_tn_marker,
        ott_tn_img,
        mark_ott_list,
        px_ott,
    ) = roc.gimmetheconf(tnconf)
    return (
        cgh_tn_marker,
        cgh_tn_img,
        tnpar,
        mark_cgh_list,
        f0,
        f1, 
        ott_tn_marker,
        ott_tn_img,
        mark_ott_list,
        px_ott,
    )


def copyConf(tnconf, tnpar):
    fromf = (
        read_ottcalib_conf.basepath + "/" + read_ottcalib_conf.fold + tnconf + ".ini"
    )
    tof = foldname.PARABOLA_REMAPPED_FOLDER + "/" + tnpar + "/" + tnconf + ".ini"
    shutil.copyfile(fromf, tof)



def view_calibration(imgott, imgpar, imgres, vm=50e-9, crpar=None, nopsd=0):
    """
    function to visualize the subtraction results.
    inputs: imgott: image of the OTT
            imgpar: registered image of the par, filtered in the case
            vm: colorbar limits
    """
    # step1: viewing the frames
    view_subplots(imgott, imgpar, imgres, crpar=None, vm=vm, nopsd=nopsd)
    view_subplots(imgott, imgpar, imgres, crpar, vm, nopsd=nopsd)


def view_subplots(imgott, imgpar, imgres, crpar=None, vm=50e-9, nopsd=0):
    if crpar is not None:
        x = crpar[0]
        y = crpar[1]
        c = crpar[2]
        oc = imgott[x : x + c, y : y + c]
        oc = zern.remov_zernike(oc, [1, 2, 3])
        pc = 2 * imgpar[x : x + c, y : y + c]
        pc = zern.removeZernike(pc, [1, 2, 3])
        rc = imgres[x : x + c, y : y + c]
        rc = zern.remove_zernike(rc, [1, 2, 3])
    else:
        oc = imgott.copy()
        pc = 2 * imgpar
        rc = imgres.copy()

    plt.figure(figsize=(18, 6))
    ax1 = plt.subplot(1, 3, 1)
    ax1.imshow(oc, vmin=-vm, vmax=vm)
    ax1.colorbar()
    ax1.set_title("OTT Image" + rmstitle(oc))

    ax2 = plt.subplot(1, 3, 2)
    ax2.imshow(pc, vmin=-vm, vmax=vm)
    ax2.colorbar()
    ax2.set_title("2Par Image" + rmstitle(pc))

    ax3 = plt.subplot(1, 3, 3)
    ax3.imshow(rc, vmin=-vm, vmax=vm)
    ax3.colorbar()
    ax3.set_title("Residue" + rmstitle(rc))

    if crpar is not None and nopsd == 0:
        norm = "ortho"
        px_ott = 0.00076
        plt.clf()
        plt.figure()
        xo, yo = th.comp_psd(oc, d=px_ott, norm=norm, verbose=True)
        xp, yp = th.comp_psd(pc, d=px_ott, norm=norm, verbose=True)
        xr, yr = th.comp_psd(rc, d=px_ott, norm=norm, verbose=True)
        plt.plot(xo[1:], yo[1:] * xo[1:], "o")
        plt.plot(xp[1:], yp[1:] * xp[1:], "x")
        plt.plot(xr[1:], yr[1:] * xr[1:], "r")
        plt.legend(["OTT", "2Par", "Res"])
        plt.xscale("log")
        plt.yscale("log")
        plt.grid()




def adjust_marker(cghf, ottf, mid, ran):
    cf = marker_remap(cghf, ottf)
    dd0 = np.sqrt((cf[0, :] - ottf[0, :]) ** 2 + (cf[1, :] - ottf[1, :]) ** 2)
    x = np.linspace(-ran, ran, 2 * ran)
    y = np.linspace(-ran, ran, 2 * ran)
    w = np.zeros((2 * ran, 2 * ran))
    for ii in np.arange(len(x)):
        for jj in np.arange(len(y)):
            tmp = cghf.copy()
            tmp[:, mid] = tmp[:, mid] + [x[ii], y[jj]]
            cf = marker_remap(tmp, ottf)
            dd = np.sqrt((cf[0, :] - ottf[0, :]) ** 2 + (cf[1, :] - ottf[1, :]) ** 2)
            w[ii, jj] = dd.std()
    a = np.array(np.where(w == np.min(w))).flatten()
    plt.imshow(w)
    plt.colorbar()
    plt.title("Marker scatter")
    print(a)
    ss = w.shape
    offs = np.array([x[a[0]], y[a[1]]])
    print(offs)
    tmp = cghf.copy()
    tmp[:, mid] = tmp[:, mid] + offs
    print("Initial pos error sigma:")
    print(dd0.std())
    print("Minimum pos error sigma:")
    print(np.min(w))
    return tmp

def rmstitle(rr):
    out = " SfE= " + str(int(rr.std() * 1e9)) + "nm"
    return out

