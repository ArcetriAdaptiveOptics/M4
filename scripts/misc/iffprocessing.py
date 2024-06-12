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

trigFrame  = 0
regFrame   = 0

class StopIteration(Exception):
    pass


def yyy(a,b)
    """
    This function ...
    Parameters
    ----------------
    Returns
    -------
    """
    return xxx

def iffRedux(flist):
    """
    Applies the differenzial algorithm...

    Parameters
    ----------


    Returns
    -------

    """



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

    return trigFrame, regFrames
