'''
this module includes the functions for the processing and analysis of the IFF data
'''
from m4.configuration.read_iffconfig import getConfig
from m4.mini_ott import timehistory as th

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




    return triggerId

