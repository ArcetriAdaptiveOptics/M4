'''
this module includes the functions for the preparation and the acquisition  IFF data
'''
from m4.configuration import read_iffconfig as rif  #contains a dictionary or similar with configuration
from m4.device import deformable_mirror as dm
from m4.iffutils import iff_acquisition_preparation as ifa
m4u = dm.M4AU()
ifc=ifa.IFFCapturePreparation(m4u)

#Template
def yyy(a,b):
    """
    This function ...
    Parameters
    ----------------
    Returns
    -------
    """
    return 


def iffDataCollection(modesList, amplitude, template=None, modalBase = None,shuffle = False):
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
    ifc._modalBaseId = modalbase
    tmh = ifc.createTimedCmdMatrixHistory(modesList,amplitude, template, shuffle )
    tn = Timestamp.now()
    _saveMatrix(filename, ifc._cmdMatrix)
    _saveMatrix(filenamex,amplitude)
    #tn = prepareIFFcollection(modesList, amplitude, nRepetitions, modalBase = None)
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


