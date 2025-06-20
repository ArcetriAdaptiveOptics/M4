"""
Author(s)
    - Pietro Ferraiuolo : written in 2024
    - Runa Briguglio : written in 2024

Description
-----------
This module contains the class that defines the M4 Adaptive Unit (M4AU) device.
"""
import os
import time
import numpy as np
from astropy.io import fits as pyfits
from m4.configuration import config_folder_names as fn
from m4.devices.base_deformable_mirror import BaseDeformableMirror
from m4.ground import logger_set_up as lsu, timestamp, read_data as rd

#put the imports from Mic Library
from Microgate.adopt.AOClient import AO_CLIENT
#this import reads automatically the HW configuration of the system, as a TN conf passed. in the TN all the relevant data are stored

#here we define the M4 class. 
#implement it properly!!!
_ts = timestamp.Timestamp()
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

    def __init__(self, tracknum=None):
        """The constructor"""
        """
        print(f"Initializing the M4AU with configuration: '{os.path.join(fn.MIRROR_FOLDER,tracknum)}'")
        self.dmConf      = os.path.join(fn.MIRROR_FOLDER,tracknum)
        self.nActs       = self._initNActuators()
        self.mirrorModes = self._initMirrorModes()
        self.actCoord    = self._initActCoord()
        self.workingActs = self._initWorkingActs()
        """
        self._aoClient = AO_CLIENT(tracknum)
        ffm = (self._aoClient.aoSystem.sysConf.gen.FFWDSvdMatrix)[0]# 
        ff  = self._aoClient.aoSystem.sysConf.gen.FFWDMatrix

        print('init the DM with no configurations')

    def getCounter(self):
        """
        Function which returns the current shape of the mirror.

        Returns
        -------
        shape: numpy.ndarray
            Current shape of the mirror.
        """
        fc = self._aoClient.getCounters()
        skipByCommand = fc.skipByCommand
        #.....
        return skipByCommand

    def get_shape(self):
        """
        Function which returns the current shape of the mirror.

        Returns
        -------
        shape: numpy.ndarray
            Current shape of the mirror.
        """
        #micLibrary.get_position()
        pos = self._aoClient.getPosition()
        return pos


    def get_force(self):
        """
        Function which returns the current force applied to the mirror.

        Returns
        -------
        force: numpy.ndarray
            Current force applied to the mirror actuators.

        """
        #micLibrary.getForce()
        force = self._aoClient.getForce()
        return force

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
        self._aoClient.mirrorCommand(cmd)
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

class AlpaoDm(BaseDeformableMirror):
    """
    Alpao interface with M4 software.
    """

    def __init__(self, ip:str, port:int):
        """The Contructor"""
        import plico_dm
        print("Ricorda di spostare le _dmCoords!")
        self._dmCoords      = {
            'dm97' : [5, 7, 9, 11],
            'dm277': [7, 9, 11, 13, 15, 17, 19],
            'dm468': [8, 12, 16, 18, 20, 20, 22, 22, 24],
            'dm820': [10, 14, 18, 20, 22, 24, 26, 28, 28, 30, 30, 32],
        }
        self._dm            = plico_dm.deformableMirror(ip, port)
        self.nActs          = self._initNactuators()
        self.mirrorModes    = None
        self.actCoord       = self._initActCoord()
        self.cmdHistory     = None
        self.baseDataPath   = fn.OPD_IMAGES_ROOT_FOLDER
        self.refAct         = 425

    def get_shape(self):
        shape = self._dm.get_shape()
        return shape
    
    def set_shape(self, cmd, differential:bool=False):
        if differential:
            shape = self._dm.get_shape()
            cmd = cmd + shape
        self._checkCmdIntegrity(cmd)
        self._dm.set_shape(cmd)

    def uploadCmdHistory(self, cmdhist):
        self.cmdHistory = cmdhist

    def runCmdHistory(self, interf=None, delay=0.2, save:str=None, differential:bool=True):
        if self.cmdHistory is None:
            raise ValueError("No Command History to run!")
        else:
            tn = _ts.now() if save is None else save
            print(f"{tn} - {self.cmdHistory.shape[-1]} images to go.")
            datafold = os.path.join(self.baseDataPath, tn)
            s = self.get_shape()
            if not os.path.exists(datafold) and interf is not None:
                os.mkdir(datafold)
            for i,cmd in enumerate(self.cmdHistory.T):
                print(f"{i+1}/{self.cmdHistory.shape[-1]}", end="\r", flush=True)
                if differential:
                    cmd = cmd+s
                self.set_shape(cmd)
                if interf is not None:
                    time.sleep(delay)
                    img = interf.acquire_phasemap()
                    path = os.path.join(datafold, f"image_{i:05d}.fits")
                    rd.save_phasemap(path, img)
        self.set_shape(s)
        return tn

    def setZeros2Acts(self):
        zero = np.zeros(self.nActs)
        self.set_shape(zero)

    def nActuators(self):
        return self.nActs
    
    def _checkCmdIntegrity(self, cmd):
        mcmd = np.max(cmd)
        if mcmd > 0.9:
            raise ValueError(f"Command value {mcmd} is greater than 1.")
        mcmd = np.min(cmd)
        if mcmd < -0.9:
            raise ValueError(f"Command value {mcmd} is smaller than -1.")
        scmd = np.std(cmd)
        if scmd > 0.5:
            raise ValueError(f"Command standard deviation {scmd} is greater than 0.1.")

    def _initNactuators(self):
        return self._dm.get_number_of_actuators()

    def _initActCoord(self):
        nacts_row_sequence = self._dmCoords[f'dm{self.nActs}']
        n_dim = nacts_row_sequence[-1]
        upper_rows = nacts_row_sequence[:-1]
        lower_rows = [l for l in reversed(upper_rows)]
        center_rows = [n_dim]*upper_rows[0]
        rows_number_of_acts = upper_rows + center_rows + lower_rows
        N_acts = sum(rows_number_of_acts)
        n_rows = len(rows_number_of_acts)
        cx = np.array([], dtype=int)
        cy = np.array([], dtype=int)
        for i in range(n_rows):
            cx = np.concatenate((cx, np.arange(rows_number_of_acts[i]) + (n_dim - rows_number_of_acts[i]) // 2))
            cy = np.concatenate((cy, np.full(rows_number_of_acts[i], i)))
        self.actCoord = np.array([cx, cy])
        return self.actCoord
