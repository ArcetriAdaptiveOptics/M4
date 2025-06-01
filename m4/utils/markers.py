import numpy as np
import os
from m4.ground import read_data
from m4.ground.timestamp import Timestamp
from m4.mini_OTT import timehistory as th
from m4.ground import geo
from m4.configuration import config_folder_names as foldname
import configparser

config = configparser.ConfigParser()
from m4.utils import image_registration_lib as imgreg
from m4.utils.parabola_identification import ParabolaActivities

pa = ParabolaActivities()

phasemapname = "surfMap.fits"
fringesname = "fringesImage.fits"
markcentername = "markerCenter.ini"
markersPath = os.path.join(foldname.BASE_PATH, foldname.MARKERS_ROOT_FOLDER)
# modificare interferometer per aggiungere loadConfiguration
markersConfig = "D:/config/20240608_negativeMarkersMask50mm.ini"  # negative mask (pass inside marker area)
markersDiam = 10
mlist0 = [2, 3, 4]
mlist1 = [
    0,
    1,
    5,
    6,
    8,
]  # removed the act out of symmetry, or the fitEllipse will find an ellipse


def acquireMarkersData(interf):
    """
    Parameters
    ----------
    interf: interferometer object
        measurement file path
    Returns
    ----------
    tn: string
        tracking number where data are saved
    """
    tn = Timestamp.now()
    fold = os.path.join(markersPath, tn)
    interf.loadConfiguration(
        markersConfig
    )  # such markersConfig is negative, i.e. mask outside the markers
    q = (
        interf.acquire_phasemap()
    )  # this is requested since the detector image has no processing to add the mask
    q = interf.intoFullFrame(q)
    ima = interf.acquire_detector()
    imamask = np.ones(ima.shape)
    imamask[ima < 0.5 * np.nanmean(ima)] = 0  # was 2x
    ima = np.ma.masked_array(ima, -1 * imamask + 1)
    # ima = np.ma.masked_array(ima, q.mask) #to be deleted
    ima = interf.intoFullFrame(ima)
    os.mkdir(fold)
    interf.save_phasemap(fold, phasemapname, q)
    interf.save_phasemap(fold, fringesname, ima)
    print(tn)
    return tn


def getParCenter(tn=None):
    """
    This function read the file with marker pattern center
    Parameters
    ----------
    tn: string
        folder where data are saved
    Returns
    ----------
    parCenter: array
        parabola center (center of markers pattern)
    """
    if tn is None:
        z = sorted(os.listdir(markersPath))[-1]
        tn = z
        print("Using the last Tracknum: " + tn)

    filename = os.path.join(markersPath, tn, markcentername)
    config.read(filename)
    cc = config["PARABOLA"]
    parCenter = np.array(cc["center"])
    return parCenter


def loadMarkersData(tn):
    """
    Parameters
    ----------
    tn: string
        folder where data are saved
    Returns
    ----------
    img: masked array
        detector image
    """
    imgname = os.path.join(markersPath, tn, fringesname)
    img = read_data.read_phasemap(imgname)  # img     = read_data.readFits_data(imgname)
    return img


def findMarkers(img):
    """
    This function find the markers in a detector image. The frame is supposed to be masked with a markers mask (individual mask diameters larger than actual markers size)
    Parameters
    ----------
    img: masked array
        masked detector image
    Returns
    ----------
    mpos: array (2xn)
        markers position (same coordinate order as the frame
    """
    mpos = pa.rawMarkersPos(img)
    npix = 3.14 * (markersDiam / 2) ** 2
    athr = 0.7
    pos = pa.filterMarkersPos(mpos, (-athr) * npix, (1 + athr) * npix)
    mpos = np.array([pos[1, :], pos[0, :]])
    return mpos


def findMarkersOrigin(pos):
    """
    This function finds the center of the markers pattern to allow identify the PAR position in the frame
    Parameters
    ----------
    pos: array (2xn)
        position of markers in the frame (same coordinate order as image)
    Returns
    ----------
    parC0: array (2x1)
        center of par
    """

    c0, axs0, r0 = pa._fitEllipse(
        pos[0, mlist1], pos[1, mlist1]
    )  # we use just the outer ring and swap coord
    c0 = np.real(c0)
    return c0


def measureMarkerPos(tn=None, interf=None):
    if tn is None:
        tn = acquireMarkersData(interf)
    print(tn)
    img = loadMarkersData(tn)
    pos = findMarkers(img)
    c0 = findMarkersOrigin(pos)
    fname = os.path.join(markersPath, tn, markcentername)
    fl = open(fname, "w")
    fl.write("[PARABOLA]" + "\n")
    data = "center  = [{:.2f},{:.2f}]".format(c0[0], c0[1])
    print(data)
    fl.write(data)
    fl.close()
    return c0


'''
def maskDetectorImage(img):
    """
    This function masks the detector image where the intensity is below a given threshold
    Parameters
    ----------
    img: masked array
        detector image
    Returns
    ----------
    imgout: masked array
        detector image masked where intensity is below threshold
    """
    mask = -1*img.mask+1
    imamask = np.ones(q.shape)
    imamask[ima <2*ima.mean()]=0
    ima = np.ma.masked_array(ima, -1*imamask+1)
    imgout = np.ma.masked_array(img, newmask)
    return imgout
'''

# how to find the M4 segment position in frame:
# Sc: segment center, coordinate of the central actuator in segment (defined)
# Scf: nominal position in the frame of the Sc posistion (this is defined by moving the truss to aim a given segment
# procedure: we acquire a number of actuators and measure their position in the frame (MAct);
# then we build the system with nominal coordinate (PAct) and measured ones (MAct).
# the current Scf_i is found with:
#   Scf_i = imgreg.marker_general_remap(MAct, PAct, Scf) (check the inputs)
