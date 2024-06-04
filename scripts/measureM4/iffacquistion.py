'''
this module includes the functions for the preparation and the acquisition  IFF data
'''
from m4.configuration import iffConfig as iffc  #contains a dictionary or similar with configuration

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


def iffCapture(first, last, amplitude, nRepetitions, modalBase = None):
    """
    This is the user-level function for the sampling of the IFF data
    Parameters
    ----------------
    first, last: int
        first and last mode index to be measured, relative to the command matrix to be used
    amplitude: float
        command amplitude
    nRepetitions: int
        number of push and pull to be collected
    modalBase: string
        identifier of the modal base to be used
    Returns
    -------
    tnif: string
        tracking number of the dataset acquired
    """
    #qui decidere se questo crea i comandi o solo fa l'acqusizione
    return tnif

def createRegistrationCommands():


