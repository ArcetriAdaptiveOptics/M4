'''
Authors
    - P. Ferriauolo
    - R. Briguglio

    Written in june 2024

------------------------

Module which contains functions for data processing and analysis of acquired M4's IFF.
'''
import os
from m4.configuration import read_iffconfig
from m4.mini_ott import timehistory as th
from m4.ground import read_data as rd
from m4.configuration import config_folder_names as foldname

imgFold = foldname.OPD_IMAGES_ROOT_FOLDER
ifFold  = foldname.IFFUNCTION_ROOT_FOLDER
# intMatFold = foldname.INTMAT_ROOT_FOLDER --- To be defined

def process(tn):
    """
    """
    regFrames, imgList = getTriggerAndRegistrationFrames(tn)
    iffRedux(imgList)
    registerRedux(regFrames)

def _ampReorganization():
    pass

def _indexReorganization():
    pass

def iffRedux(flist):
    """
    Apply the differential algorythm.

    Parameters
    ----------
    flist : str / ArrayLike
        Images file list or folder tracking number

    Returns
    -------
    cube ?

    """
    _,_,_,template  = read_iffconfig('IFFUNC')
    nPushPull = len(template)
    indexList = _indexReorganization()    # pass
    amplitude = _ampReorganization()      # pass

    for i in range(nModes):
        print(i)
        for k in range(nPushPull)):
            p = nPushPull * i + k
            n = where[p]
            mis_amp = k * indexingList.shape[1] + n

        image = # caricata con la maschera da flist

# Algorimo centrale
        for p in range(1, len(flist)):
            opd2add = image[p]*template[p] + image[p-1]*template[p-1]
            master_mask2add = np.ma.mask_or(image[p].mask, image[p-1].mask)
            if p==1:
                master_mask = master_mask2add
            else:
                master_mask = np.na.mask_or(master_mask, master_mask2add)
            image += opd2add

        image = np.ma.masked_array(image, mask=master_mask)
        norm_image = image / (2*amp_reorg[mis_amp] * (template.shape[1]-1))
        rd.save_phasemap(os.path.join(intMatFold, 'mode_{:5d}.fits'.format(i)), norm_image)


def getTriggerAndRegistrationFrames(flist, amplitude=None)
    """
    This function identifies the triggering frame and the registration frames in the data frame sequence.

    Parameters
    ----------
    flist : str / list
        Complete list of the files. Can be either the folder's tracking number, containing the images, or the list of images itself.
    amplitude: float
        Amplitude of the trigger command applied. If no value is passed, it will be loaded from the 'iffconfig.ini' configuration file.

    Returns
    -------
    trigFrame : int
        Index of the frame containing the trigger. The subsequent frame (triggerId + 1) is the first 'useful' frame in the time history.
    regFrames : int | ArrayLike
        List containing the first and last indices of the registration frames. The subsequent frame is the start of the IFF acquisition commanded matrix history.
    """
    if isinstance(flist, str):
        filelist = th.findtracknum(flist)

    configfile = os.path.join(os.path.dirname(filelist[0]), 'iffconfig.iini')

    if amplitude is None
        _,_,amplitude,_ = read.iffconfig.getConfig('TRIGGER', bpath=configfile)

    global trigFrame
    global regFrames

    try:
        for ii in range(len(filelist)):
            f0 = rd.readFits_maskedImage(filelist[ii])
            f1 = rd.readFits_maskedImage(filelist[ii+1])

            diff = th.removeZernike(f1 - f0)
            std_check = np.std(diff)

            if (std_check>(amp/2)): # condizione
                trigFrame = ii
                raise StopIteration
    except StopIteration:
        pass

    timing = read_iffconfig.getTiming()
    regZeros, regModes, _, regTemplate = read.iffconfig.getConfig('REGISTRATION', bpath=configfile)

    regStart  = trigFrame + regZeros*timing
    regEnd    = len(regModes)*len(regTemplate)*timing
    regFrames = [regStart, regEnd]

    imgList = filelist[(regFrames[1]+1):]

    return regFrames, imgList

class StopIteration(Exception):
    pass