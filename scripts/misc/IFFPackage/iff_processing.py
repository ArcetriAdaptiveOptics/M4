"""
Authors
    - P. Ferriauolo
    - R. Briguglio

    Written in june 2024

------------------------

Module which contains functions for data processing and analysis of M4's IFF.
"""
import os
import numpy as np
from astropy.io import fits as pyfits
from m4.ground import timestamp
from m4.configuration import read_iffconfig
from m4.mini_OTT import timehistory as th
from m4.ground import read_data as rd
from m4.configuration import config_folder_names as fn
from m4.devices import deformable_mirror as dm
from scripts.misc.IFFPackage import iff_acquisition_preparation

ifa = iff_acquisition_preparation.IFFCapturePreparation(dm.M4AU())

imgFold     = fn.OPD_IMAGES_ROOT_FOLDER
ifFold      = fn.IFFUNCTIONS_ROOT_FOLDER
intMatFold  = fn.INTMAT_ROOT_FOLDER
confFold    = fn.CONFIGURATION_ROOT_FOLDER

def process(tn):
    """
    

    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    regFrames, imgList = getTriggerAndRegistrationFrames(tn)
    iffRedux(imgList)
    registrationRedux(regFrames)

def iffRedux(tn):
    """
    Processes 4D images to produce FITS images
    
    This function retrieves the 4D images acquired and pre-processed from the
    interferometer, applies the differential algorythm, and produces a fits image
    for each measured mode.

    Parameters
    ----------
    tn : str 
        Tracking number in the OPDImages folder where the acquired data is stored

    Returns
    -------
    

    """
    _,nModes,_,template,_  = read_iffconfig.getConfig('IFFUNC')
    nPushPull = len(template)
    indexingList = _indexReorganization()    # to be implemented
    amplitude = _ampReorganization()      # to be implemented
    filelist = os.listdir(os.path.join(imgFold, tn))
    for i in range(nModes):
        print(i)
        for k in range(nPushPull):
            p = nPushPull * i + k
            n = indexingList[p]
            mis_amp = k * indexingList.shape[1] + n
            mis = mis_amp*template.shape[0]
            # Crea il pacchetto di immagini del modo 'i', contenente nPushPull images
            file_name = filelist[mis]
            image_list = []
            for l in range(0, template.shape[0]):
                file_name = filelist[mis+l]
                ima = rd.read_phasemap(file_name)
                image_list.append(ima)
            image = np.zeros((ima.shape[0], ima.shape[1]))
        # Algorimo differenziale
        for p in range(1, len(filelist)):
            opd2add = image[p]*template[p] + image[p-1]*template[p-1]
            master_mask2add = np.ma.mask_or(image[p].mask, image[p-1].mask)
            if p==1:
                master_mask = master_mask2add
            else:
                master_mask = np.na.mask_or(master_mask, master_mask2add)
            image += opd2add
        image = np.ma.masked_array(image, mask=master_mask)
        norm_image = image / (2*amplitude[mis_amp] * (template.shape[1]-1))
        fold = os.path.join(ifFold, tn)
        rd.save_phasemap(os.path.join(fold, 'mode_{:5d}.fits'.format(i)), norm_image)
        
def createCube(tn):
    """
    

    Parameters
    ----------
    tn : TYPE
        DESCRIPTION.

    Returns
    -------
    cube : TYPE
        DESCRIPTION.

    """
    filelist = th.fileList(tn)
    cube_list = []
    for imgfits in filelist:
        image = rd.readFits_maskedImage(imgfits)
        cube_list.append(image)
    cube = np.ma.dstack(cube_list)
    path = os.path.join(fn.INTMAT_ROOT_FOLDER, tn)
    save_name = os.path.join(path, 'IMCube.fits')
    rd.save_phasemap(save_name, cube, isCube=True)
    print("Cube saved in '{}'".format(path))
    return cube #!!! return?

def stackCubes(tnlist):
    """
    

    Parameters
    ----------
    tnlist : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    search_fold = fn.IFFUNCTIONS_ROOT_FOLDER
    cubelist = []
    for tn in tnlist:
        fold = os.path.join(search_fold, tn)
        cube_name = os.path.join(fold, 'IMCube.fits')
        hdul = pyfits.open(cube_name)
        cubelist.append(hdul[0].data)
        hdul.close()
    stacked_cube = np.ma.dstack(cubelist)
    new_tn = timestamp.Timestamp.now()
    save_fold = os.path.join(fn.INTMAT_ROOT_FOLDER, new_tn)
    rd.save_phasemap(os.path.join(save_fold, 'IMCube.fits'), stacked_cube, isCube=True)
    print("Stacked cube saved in {}".format(save_fold))
    return stacked_cube #!!! Return?
    

def registrationRedux(regFrames):
    """
    

    Parameters
    ----------
    regFrames : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    pass

def getTriggerAndRegistrationFrames(flist, amplitude=None):
    """
    This function identifies the triggering frame and the registration frames in
    the data frame sequence.

    Parameters
    ----------
    flist : str / list
        Complete list of the files. Can be either the folder's tracking number,
        containing the images, or the list of images itself.
    amplitude: float
        Amplitude of the trigger command applied. If no value is passed, it will
        be loaded from the 'iffconfig.ini' configuration file.

    Returns
    -------
    trigFrame : int
        Index of the frame containing the trigger. The subsequent frame (triggerId + 1) 
        is the first 'useful' frame in the time history.
    regFrames : int | ArrayLike
        List containing the first and last indices of the registration frames. 
        The subsequent frame is the start of the IFF acquisition commanded matrix
        history.
    """
    if isinstance(flist, str):
        filelist = th.fileList(flist)
    if amplitude is None:
        _,_,amplitude,_, _ = read_iffconfig.getConfig('TRIGGER')
    try:
        for ii in range(len(filelist)):
            f0 = rd.readFits_maskedImage(filelist[ii])
            f1 = rd.readFits_maskedImage(filelist[ii+1])
            diff = th.removeZernike(f1 - f0)
            std_check = np.std(diff)
            if std_check>(amplitude/2): # condizione -- check
                trigFrame = ii
                raise StopIteration
    except StopIteration:
        pass
    timing = read_iffconfig.getTiming()
    regZeros, regModes, _, regTemplate, _ = read_iffconfig.getConfig('REGISTRATION')
    regStart  = trigFrame + regZeros*timing
    regEnd    = len(regModes)*len(regTemplate)*timing
    regFrames = [regStart, regEnd]
    imgList = filelist[(regFrames[1]+1):]
    return regFrames, imgList

def _ampReorganization():
    """
    

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    indexList = ifa.getIndexingList()
    return indexList

def _indexReorganization():
    """
    

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    amps = ifa.getAmplitude()
    return amps

class StopIteration(Exception):
    pass
