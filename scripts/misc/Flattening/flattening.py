"""
Author(s) 
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
import os
import numpy as np
from m4.dmutils import iff_processing as ifp
from m4.ground import read_data as rd
from m4.analyzers import compute_reconstructor as crec

class Flattening:
    """
    This class xxxxx

    Methods
    ======
    loadIntMat
    computeRecMat
    flat
    ...
    """
    def __init__(self, tn):
        """The Constructor"""
        self._tn            = tn
        self._path          = os.path.join(ifp.intMatFold, self._tn)
        self._intMat        = self._loadIntMat()
        self._cmdMat        = self._loadCmdMat()
        self._rec           = crec.ComputeReconstructor(self._intMat)
        self._recMat        = None
        self._frameCenter   = None
        self._shape2flat    = None # Immagine acquisita o fornita, come? 
        self._flatOffset    = None
        self._cavityOffset  = None
        self._synthFlat     = None
        self._flatResidue   = None
        self._flatResult    = None
        self._flatteningModes = None

    def computeFlatCmd(self):
        """
        routine which computes the command to apply to flatten the input shape.

        Returns
        -------
        flat_cmd : ndarray
            Flat command.
        """
        # Logica
        cmd_amp = -np.dot(self._recMat.flatten(), self._shape2flat.compressed())
        flat_cmd = np.dot(self._cmdMat, cmd_amp)
        self._flatResult = flat_cmd
        return flat_cmd

    def load_image2shape(self, img, compute_rec:bool=True):
        """
        (Re)Loader for the image to flatten.

        Parameters
        ----------
        img : MaskedArray
            Image to flatten.
        compute_rec : bool, optional
            Wether to direclty compute the reconstructor with the imput image or
            not. The default is True.
        """
        self._shape2flat = img
        self._rec.loadShape2Flat(img)
        if compute_rec:
            print("Computing recontruction matrix...")
            self._recMat = self._rec.run()

    def _loadIntCube(self):
        """
        Interaction cube loader
        """
        cube = os.path.join(self._path, ifp.cubeFile)
        with open(os.path.join(self._path, ifp.flagFile), 'r', encoding='utf-8') as f:
            flag = f.read()
        if ' filtered ' in flag:
            intCube = rd.read_phasemap(cube)
        else:
            intCube = ifp.filterZernikeCube(self._tn)
        return intCube

    def _loadCmdMat(self):
        """
        Command matrix loader. It loads the saved command matrix of the loaded
        cube.

        Returns
        -------
        cmdMat : ndarray
            Command matrix of the cube, saved in the tn path.
        """
        cmdMat = rd.readFits_data(os.path.join(self._path, ifp.cmdMatFile))
        return cmdMat

    def _loadFrameCenter(self):
        """
        Center frame loader, useful for image registration.

        Returns
        -------
        frame_center : TYPE
            DESCRIPTION.

        """
        frame_center = rd.readFits_data('data')
        return frame_center

    def _registerShape(self, shape):
        xxx=None
        dp = ifp.findFrameOffset(self._tn,xxx)
        #cannot work. we should create a dedicated function, not necessarily linked to IFF or flattening
        return dp
