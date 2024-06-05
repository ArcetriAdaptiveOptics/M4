'''
this module includes the functions for the preparation and the acquisition  IFF data
'''
from m4.configuration import iffConfig as iffc  #contains a dictionary or similar with configuration
from m4.device import deformable_mirror as dm

'''Template
def yyy(a,b)
    """
    This function ...
    Parameters
    ----------------
    Returns
    -------
    """
    return xxx
'''

def iffDataCollection(modesList, amplitude, nRepetitions, modalBase = None):
    """
    This is the user-level function for the sampling of the IFF data
    Parameters
    ----------------
    modesList: int | list, array like
        list of modes index to be measured, relative to the command matrix to be used
    amplitude: float
        command amplitude
    nRepetitions: int
        number of push and pull to be collected
    modalBase: string
        identifier of the modal base to be used
    Returns
    -------
    tn: string
        tracking number of the dataset acquired
    """
    tn = prepareIFFcollection(modesList, amplitude, nRepetitions, modalBase = None)
    iffCapture(tn)
    return tn



def iffCapture(tn):
    """
    This function manages the interfacing equence for collecting the IFF data
    Parameters
    ----------------
    tn: string
        the tracking number in the xxx folder where the cmd history is saved
    Returns
    -------
    """

    cmdHist = getCmdHist(tn)
    dm.uploadCmdHist(cmdHist)
    dm.runCmdHist()
    print('Now launching the acquisition sequence')
    start4DAcq(tn)
    print('Acquisition completed. Dataset tracknum:')
    print(tn)


def getCmdHist(tn):
    """
    This function ...
    Parameters
    ----------------
    Returns
    -------
    """
    cmdHist = [1,2,3]
    return cmdHist


