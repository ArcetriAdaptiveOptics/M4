'''
Authors
  - C. Selmi: written in ?
'''
import logging
import numpy as np
from m4.devices.base_deformable_mirror import BaseDeformableMirror
#put the imports from Mic Library


def mirrorCommand(cmd, segment=None):
    """
    Function for the user-level application of commands to the mirror
    Parameters
    ----------------
    cmd: float | list, array
        command to be applied [m], differential wrt the current bias command. may be full command (5000+) or single segment command (892)
    segment: int
        id of the segment to apply a given command
    Returns
    -------
    """
    print('Command applied')

def uploadCmdHist(cmdHist, timeInfo=None):
    """
    This function ...
    Parameters
    ----------------
    Returns
    -------
    """
    print('Command history uploaded')

def runCmdHist():
    """
    This function ...
    Parameters
    ----------------
    Returns
    -------
    """
    print('Command history running...')





