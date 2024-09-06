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

    def create_interaction_matrix(self, cmdMat, zernValues, cmdAmp):
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
    def __init__(self):
        """The Constructor"""

    def apply_command(self, fullCmd):
        """


        Parameters
        ----------
        fullCmd : TYPE
            DESCRIPTION.

        """

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
