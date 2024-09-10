"""
Author(s):
    - Pietro Ferraiuolo : written in 2024

Description
-----------

How to Use it
-------------
"""
import numpy as np
from scripts.ottalign import _m4ac as mac # to change
from m4.ground import zernike as zern, read_data as rd

class AlignmentCorrection():
    """
    """
    def __init__(self, ott):
        """The Constructor"""
        self.ott        = ott
        self.cmdMat     = rd.read_phasemap('/home/pietrof/git/M4/scripts/ottalign/cmdMat.fits') #temporary
        self._devCalls  = mac.dev_calls
        self._devName   = mac.names
        self._devFncs   = self._get_devices_callables()
        self._dof       = mac.dof
        self._dofTot    = mac.cmdDof
        self._idx       = mac.slices

    def apply_command(self, fullCmd):
        """
        Function which takes 

        Parameters
        ----------
        fullCmd : list or ndarray
            Full command of the interaction matrix which command all device's 
            available motors..
        """
        device_commands = self._extract_cmd(fullCmd)
        for cmd,fnc,dev in zip(device_commands,self._devFncs,self._devName):
            if np.sum(cmd)!=0.:
                try:
                    print(f"Commanding {cmd} to {dev}")
                    fnc(cmd)
                except Exception as e:
                    print(f"Qualcosa non va con {dev}:\n", e)
            else:
                print(f'skipping null command for {dev}')

    def _extract_cmd(self, fullCmd):
        """
        Function to extract the relevant input command for each device from the
        full command of the interaction matrix, creating a readeble command for
        each device.

        Parameters
        ----------
        fullCmd : list or ndarray
            Full command of the interaction matrix which command all device's 
            available motors.
        """
        commands = []
        for d in range(len(self._dof)):
            dev_cmd = np.zeros(self._dofTot)
            dev_idx = fullCmd[self._idx[d]]
            for i,idx in enumerate(dev_idx):
                dev_cmd[self._dof[d][i]] = idx
            commands.append(dev_cmd)
        return commands

    def _get_devices_callables(self):
        """
        Returns a list of callables for the instanced object, taken from the 
        configuration .py file

        Returns
        -------
        functions : list
            List of callables, which interacts with the input object of the class.
        """
        functions = []
        for dev_call in self._devCalls:
            obj, *methods = dev_call.split('.')
            call = getattr(self.ott, obj)
            for method in methods:
                call = getattr(call, method)
            functions.append(call)
        return functions

class AlignmentCalibration():
    
    def __init__(self, acquisition_device):
        """the constructor"""
        self.ccd        = acquisition_device
        self._ccdCalls  = mac.ccd_calls
        self._ccdFncs   = self._get_ccd_callables()
        
    def image_acquisition(self, cmdAmp, template:list=None, n_frames=None):
        """
        

        Parameters
        ----------
        cmdAmp : TYPE
            DESCRIPTION.
        template : list, optional
            DESCRIPTION. The default is None.
        n_frames : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """

    def _get_ccd_callables(self):
        """
        Returns a list of callables for the instanced object, taken from the 
        configuration .py file
    
        Returns
        -------
        functions : list
            List of callables, which interacts with the input object of the class.
        """
        functions = []
        for dev_call in self._ccdCalls:
            obj, *methods = dev_call.split('.')
            call = getattr(self.ccd, obj)
            for method in methods:
                call = getattr(call, method)
            functions.append(call)
        return functions
