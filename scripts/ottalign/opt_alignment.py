"""
Author(s):
    - Pietro Ferraiuolo : written in 2024

Description
-----------

How to Use it
-------------
"""
import os, numpy as np
from m4.configuration.ott_parameters import OttParameters
from m4.ground import zernike as zern, read_data as rd
op = OttParameters()
dof = {
       'Parabola': op.PARABOLA_DOF,
       'Rm': op.RM_DOF,
       'M4': op.M4_DOF
       }

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

class AlignmentCorrection():
    """
    """
    def __init__(self, ott):
        """The Constructor"""
        self.ott = ott

    def apply_command(self, fullCmd):
        """


        Parameters
        ----------
        fullCmd : TYPE
            DESCRIPTION.

        """
        # iterare per tutti i comandi della matrice? No, non ha senso (?)
        par_cmd = self._get_par_cmd(fullCmd)
        rm_cmd  = self._get_rm_cmd(fullCmd)
        m4_cmd  = self._get_m4_cmd(fullCmd)
        cmd = [par_cmd, rm_cmd, m4_cmd]
        dev = [self.ott.parabola, self.ott.referenceMirror, self.ott.dm] #provvisional
        for command,device in zip(cmd,dev):
            if np.sum(device) != 0.:
                device.setPosition(command)
            else:
                print(f"Null command for {device}: skipping...")
        text = f"""
Parabola position : {dev[0].getPosition()}
Reference Mirror position : {dev[1].getPosition()}
M4 position : {dev[2].getPosition()}
"""
        print(text)

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

    def _get_par_cmd(self, cmd_vect):
        """
        Extract the parabola command out the full command vector.

        Parameters
        ----------
        cmd_vect : ndarray
            Full command vector.
        """
        rm = cmd_vect[3:5]
        rm_cmd = np.zeros(6)
        for i,idx in enumerate(rm):
            rm_cmd[dof['Rm'][i]] = idx
        return rm_cmd

    def _get_rm_cmd(self, cmd_vect):
        """
        Extract the parabola command out the full command vector.

        Parameters
        ----------
        cmd_vect : ndarray
            Full command vector.
        """
        m4 = cmd_vect[5:]
        m4_cmd = np.zeros(6)
        for i,idx in enumerate(m4):
            m4_cmd[dof['M4'][i]] = idx
        return m4_cmd

    def _get_m4_cmd(self, cmd_vect):
        """
        Extract the parabola command out the full command vector.

        Parameters
        ----------
        cmd_vect : ndarray
            Full command vector.
        """
        par = cmd_vect[:3]
        par_cmd = np.zeros(6)
        for i,idx in enumerate(par):
            par_cmd[dof['Parabola'][i]] = idx
        return par_cmd

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
