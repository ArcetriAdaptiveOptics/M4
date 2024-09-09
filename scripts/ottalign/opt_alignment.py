"""
Author(s):
    - Pietro Ferraiuolo : written in 2024

Description
-----------

How to Use it
-------------
"""
import os
import numpy as np
from script.ottalign import m4_alignment_configuration as mac # to change
from m4.ground import zernike as zern, read_data as rd

class AlignmentCorrection():
    """
    """
    def __init__(self):
        """The Constructor"""
        self.ott        = [mac.ottpar, mac.ottrm, mac.ottm4]
        self.dof        = mac.dof
        self.dof_tot    = mac.cmdDof
        self.idx        = mac.slices

    def apply_command(self, fullCmd):
        """


        Parameters
        ----------
        fullCmd : TYPE
            DESCRIPTION.

        """
        device_commands = self._extract_cmd(fullCmd)
        for cmd,dev in zip(device_commands,self.ott):
            if np.sum(cmd)!=0.:
                dev(cmd)


    def command_creation(self, redCmdVect, dof):
        """


        Parameters
        ----------
        redCmdVect : TYPE
            DESCRIPTION.
        dof : TYPE
            DESCRIPTION.

        Returns
        -------
        fullCmd
            DESCRIPTION.

        """

    def _extract_cmd(self, fullCmd):
        """


        Parameters
        ----------
        fullCmd : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        commands = []
        for d in range(len(self.dof)):
            dev_cmd = np.zeros(self.dof_tot)
            dev_idx = fullCmd[self.idx[d]]
            for i,idx in enumerate(dev_idx):
                dev_cmd[self.dof[d][i]] = idx
            commands.append(dev_cmd)
        return commands

    def _get_reduced_cmd_vector(self, zernVect, recMat):
        """


        Parameters
        ----------
        zernVect : TYPE
            DESCRIPTION.
        recMat : TYPE
            DESCRIPTION.

        Returns
        -------
        redCmdVect
            DESCRIPTION.

        """

    def _images_acquisition(self, n_frames):
        """


        Parameters
        ----------
        n_frames : TYPE
            DESCRIPTION.

        Returns
        -------
        images
            DESCRIPTION.

        """

    def _load_data(self, tracknum):
        """


        Parameters
        ----------
        tracknum : TYPE
            DESCRIPTION.

        """

    def _zernike_routine(self, images, zernikeList, fitArea):
        """


        Parameters
        ----------
        images : TYPE
            DESCRIPTION.
        zernikeList : TYPE
            DESCRIPTION.
        fitArea : TYPE
            DESCRIPTION.

        Returns
        -------
        zernValues
            DESCRIPTION.

        """

    def intmat_cropping(self, intmat, cmdMat, zernList, dof):
        """


        Parameters
        ----------
        intmat : TYPE
            DESCRIPTION.
        cmdMat : TYPE
            DESCRIPTION.
        zernList : TYPE
            DESCRIPTION.
        dof : TYPE
            DESCRIPTION.

        Returns
        -------
        intmat
            DESCRIPTION.

        """

    def _matrix_inversion(self, intmat):
        """


        Parameters
        ----------
        intmat : TYPE
            DESCRIPTION.

        Returns
        -------
        recMat
            DESCRIPTION.

        """


class AlignmentCalibration():
    """
    """
    def __init__(self):
        """The Constructor"""

    def image_acquisition(self, cmdMat, template, n_frames, cmdAmp):
        """


        Parameters
        ----------
        cmdMat : TYPE
            DESCRIPTION.
        template : TYPE
            DESCRIPTION.
        n_frames : TYPE
            DESCRIPTION.
        cmdAmp : TYPE
            DESCRIPTION.

        Returns
        -------
        images
            DESCRIPTION.

        """

    def create_Cmatrix(self, cmdMat, zernValues, cmdAmp):
        """


        Parameters
        ----------
        cmdMat : TYPE
            DESCRIPTION.
        zernValues : TYPE
            DESCRIPTION.
        cmdAmp : TYPE
            DESCRIPTION.

        Returns
        -------
        intMat
            DESCRIPTION.

        """

    def _zernike_routine(self, images, zernikeList, fitArea):
        """


        Parameters
        ----------
        images : TYPE
            DESCRIPTION.
        zernikeList : TYPE
            DESCRIPTION.
        fitArea : TYPE
            DESCRIPTION.

        Returns
        -------
        zernValues
            DESCRIPTION.

        """

    def _save_data(self, data):
        """


        Parameters
        ----------
        data : TYPE
            DESCRIPTION.

        """
