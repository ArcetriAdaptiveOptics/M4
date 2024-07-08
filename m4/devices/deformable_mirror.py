'''
Authors
  - C. Selmi: written in ?
'''
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

class AdOpticaDM(BaseDeformableMirror):  #ereditare BaseDeformableMirror

    def __init__(self,tracknum):
        print('Initializing the M4AU with configuration: '+fn.DM_CONFIGURATION_ID)
        self.dmConf           = os.path.join(fn.MIRROR_FOLDER,tracknum)
        self.nActs       = self._initNActuators()
        self.mirrorModes = self._initMirrorModes()
        self.actCoord    = self._initActCoord()
        self.workingActs = self._initWorkingActs()

    def nActuators(self):
        return self.nActs


    def _initNActuators(self):
        fname = os.path.join(self.dmConf,nActFile)
        f = open(fname,'r')
        nact = int(f.read())
        f.close
        return nact

    def _initMirrorModes(self):
        '''
        Creation of the mirror modes variable

        '''
        fname = os.path.join(self.dmConf,mirrorModesFile)
        if os.path.exists(fname):
            print('Initializing mirror modes from data: nact x nmodes')
            hdu = pyfits.open( fname)
            mirrorModes = hdu[0].data
        else:
            print('Initializing analytical modal base (identity, or zonal matrix')
            mirrorModes = np.eye(self.nActs)
        #nActs = np.shape(cmdMat)[0]
        return mirrorModes

    def _initWorkingActs(self):
        '''
        Reading the list of working actuators
        '''
        fname = os.path.join(self.dmConf,mirrorModesFile)
        if os.path.exists(fname):
            print('Initializing mirror modes from data: nact x nmodes')
            hdu = pyfits.open( fname)
            mirrorModes = hdu[0].data
        else:
            print('Initializing analytical modal base (identity, or zonal matrix')
            mirrorModes = np.eye(self.nActs)
        #nActs = np.shape(cmdMat)[0]
        return mirrorModes
    
    def _initActCoord(self):
        '''
        Reading the actuators coordinate from file
        '''
        fname = os.path.join(self.dmConf,actCoordFile)
        hdu = pyfits.open(fname)
        actCoord = hdu[0].data
        #nActs = np.shape(cmdMat)[0]
        return actCoord
    #def _mirrorConfFolder(self):
        #basef = fn.MIRROR_FOLDER
        #conffolder = os.path.join(basef,tn)
        #return conffolder

    def get_shape(self):
        #micLibrary.get_position()
        shape = 3
        return shape
        

    def get_force(self):
        '''
        '''
        #micLibrary.getForce()
        force = 3
        return force


    def set_shape(self, cmd):#cmd, segment=None):
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
        #micLibrary.mirrorCommand(cmd)
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





