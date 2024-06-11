'''
this module includes the functions for the processing and analysis of the IFF data
'''
from m4.configuration.read_iffconfig import getConfig
from m4.mini_ott import timehistory as th

nTrigFrame  = 0
nRegFrame   = 0

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

def ifRedux(flist,

def getRegistrationFrames(flist, triggerId)
    """
    This function ...
    Parameters
    ----------------
    Returns
    -------
    """

    return xxx

def findFirstFrame(flist, amplitude) # flist -> tn? / amplitude = iffconfig.ini?
    """
    This function identifies the triggering frame in the frames sequence

    Parameters
    ----------
    flist : 
        Complete list of the files
    amplitude: float
        Amplitude of the trigger command applied

    Returns
    -------
    triggerId: int
        Index of the frame containing the trigger. The subsequent frame (triggerId + 1) is the first 'useful' frame in the time history.
    """
    filelist = th.findtracknum(flist)
    _,_,amplitude,_ = getConfig(confile, 'TRIGGER') # poi si ragiona sul come il confile

    try:
        for file in filelist:
            with pyfits.open(file) as hdu:
                image = hdu[0].data
                mask  = hdu[1].data.astype(bool)
                count += 1
                if :#condizione
                    thefile = file
                    trigger = np.ma.masked_array(image, mask)
                    raise StopIteration
    except StopIteration:
        pass

    global nTrigFrame
    nTrigFrame = count

    return triggerId
