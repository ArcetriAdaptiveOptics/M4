"""
Author(s)
    - Pietro Ferraiuolo : written in 2024
    - Runa Briguglio : written in 2024

Description
-----------
This module contains the class that defines the M4 Adaptive Unit (M4AU) device.
"""
import os
import logging
import numpy as np
from m4.devices.base_deformable_mirror import BaseDeformableMirror
from m4.configuration import config_folder_names as fn
from astropy.io import fits as pyfits

#put the imports from Mic Library

#here we define the M4 class. 
#implement it properly!!!

mirrorModesFile = 'ff_v_matrix.fits'
ffFile          = 'ff_matrix.fits'
actCoordFile    = 'ActuatorCoordinates.fits'
nActFile        = 'nActuators.dat'

class AdOpticaDM(BaseDeformableMirror):  #ereditare BaseDeformableMirror  #forse passare un TN come configurazione del DM??
    """
    Class that defines the M4 Adaptive Unit (M4AU) device.
    
    Attributes
    ----------
    dmConf: str
        Configuration folder of the M4AU.
    nActs: int
        Number of actuators of the M4AU.
    mirrorModes: numpy.ndarray
        Mirror Modes Matrix.
    actCoord: numpy.ndarray
        Actuator Coordinates Matrix.
    workingActs: numpy.ndarray
        Working Actuators Matrix.
        
    Methods
    -------
    set_shape(cmd, segment=None)
        User-level function for the application of commands to the mirror.
    uploadCmdHist(cmdHist, timeInfo=None)
        This function loads the command matrix history into the DM.
    runCmdHist()
        This function makes the dm run the command matrix history.
    get_shape()
        Function which returns the current shape of the mirror.
    get_force()
        Function which returns the current force applied to the mirror.
    
    """

    def __init__(self, tracknum):
        """The constructor"""
        print(f"Initializing the M4AU with configuration: '{os.path.join(fn.MIRROR_FOLDER,tracknum)}'")
        self.dmConf      = os.path.join(fn.MIRROR_FOLDER,tracknum)
        self.nActs       = self._initNActuators()
        self.mirrorModes = self._initMirrorModes()
        self.actCoord    = self._initActCoord()
        self.workingActs = self._initWorkingActs()


    def set_shape(self, cmd):#cmd, segment=None):
        """
        User-level function for the application of commands to the mirror.

        Parameters
        ----------------
        cmd: float | list, array
            command to be applied [m], differential wrt the current bias 
            command. may be full command (5000+) or single segment command (892)
        segment: int
            Id of the segment to which apply the given command.

        Returns
        -------

        """
        #micLibrary.mirrorCommand(cmd)
        print('Command applied')


    def uploadCmdHist(cmdHist, timeInfo=None):
        """
        This function loads the command matrix history into the DM.

        Parameters
        ----------
        cmdHist: numpy.ndarray
            Command Matrix History.
        timeInfo: int
            Timing information of the command matrix history.

        Returns
        -------

        """
        print('Command history uploaded')


    def runCmdHist():
        """
        This function makes the dm run the command matrix history.

        Parameters
        ----------

        Returns
        -------

        """
        print('Command history running...')


    def get_shape(self):
        """
        Function which returns the current shape of the mirror.
        
        Returns
        -------
        shape: numpy.ndarray
            Current shape of the mirror.
        """
        #micLibrary.get_position()
        shape = 3
        return shape
        

    def get_force(self):
        """
        Function which returns the current force applied to the mirror.

        Returns
        -------
        force: numpy.ndarray
            Current force applied to the mirror actuators.

        """
        #micLibrary.getForce()
        force = 3
        return force


    def _initNActuators(self):
        """
        Function which reads the number of actuators of the DM from a configuration
        file.

        Returns
        -------
        nact: int
            number of actuators of the DM.
        """
        fname = open(os.path.join(self.dmConf, nActFile),'r')
        with open(fname,'r') as f:
            nact = int(f.read())
        return nact


    def _initMirrorModes(self):
        """
        Function which initialize the mirror modes by reading from a fits file.

        Returns
        -------
        mirrorModes: numpy.ndarray
            Mirror Modes Matrix.
        """
        fname = os.path.join(self.dmConf, mirrorModesFile)
        if os.path.exists(fname):
            print('Initializing Mirror Modes')
            with pyfits.open(fname) as hdu:
                mirrorModes = hdu[0].data
        else:
            print('Initializing Analytical Modal Base')
            mirrorModes = np.eye(self.nActs)
        #nActs = np.shape(cmdMat)[0]
        return mirrorModes


    def _initWorkingActs(self):
        """
        Function which initialize the working actuators by reading
        a list from a fits file.

        Returns
        -------
        workingActs: numpy.ndarray
            Working Actuators Matrix.
        """
        fname = os.path.join(self.dmConf, mirrorModesFile)
        if os.path.exists(fname):
            with pyfits.open(fname) as hdu:
                workingActs = hdu[0].data
        else:
            workingActs = np.eye(self.nActs)
        return workingActs
    

    def _initActCoord(self):
        '''
        Reading the actuators coordinate from file
        '''
        fname = os.path.join(self.dmConf, actCoordFile)
        with pyfits.open(fname) as hdu:
            actCoord = hdu[0].data
        return actCoord
    
    #def _mirrorConfFolder(self):
        #basef = fn.MIRROR_FOLDER
        #conffolder = os.path.join(basef,tn)
        #return conffolder
