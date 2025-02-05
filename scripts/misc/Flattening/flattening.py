"""
Author(s) 
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
Module containing the class which computes the flattening command for a deformable
mirror, given an imput shape and a (filtered) interaction cube.

From the loaded tracking number (tn) the interaction cube will be loaded (and
filtered, if it's not already) from which the interaction matrix will be computed.
If an image to shape is provided on class instance, then the reconstructor will
be automatically computed, while if not, the load_img2shape methos is available
to upload a shape from which compute the reconstructor.

How to Use it
=============
Instancing the class only with the tn of the interaction cube

    >>> from m4.__ import flattening as flt
    >>> tn = '20240906_110000' # example tn
    >>> flat = flt.Flattening(tn)
    >>> # say we have acquired an image
    >>> img = interf.acquire_phasemap()
    >>> flat.load_image2shape(img)
    'Computing reconstruction matrix...'

all is ready to compute the flat command, by simply running the method

    >>> flatCmd = flat.computeFlatCmd()
"""
import os
import numpy as np
from m4.ground import read_data as rd
from m4.dmutils import iff_processing as ifp
from m4.analyzers import compute_reconstructor as crec

class Flattening:
    """
    Class which handles the flattening command computation
    
    Public Methods
    -------
    computeFlatCmd : 
        Method which computes the flattening command to apply to a given shape,
        which must be already in memory, through the class instancing or the 
        load_img2shape method
    
    load_image2shape : 
        method to (re)upload and image to shape in the class, after which the
        reconstructor will be automatically computed for it.
    """
    def __init__(self, tn, img2flatten=None):
        """The Constructor"""
        self.flatCmd        = None
        self.shape2flat     = self.load_image2shape(img2flatten) \
                                           if img2flatten is not None else None
        self._tn            = tn
        self._path          = os.path.join(ifp.intMatFold, self._tn)
        self._intCube       = self._loadIntCube()
        self._cmdMat        = self._loadCmdMat()
        self._rec           = self._loadReconstructor(self._intCube)
        self._recMat        = None
        self._frameCenter   = None
        self._flatOffset    = None
        self._cavityOffset  = None
        self._synthFlat     = None
        self._flatResidue   = None
        self._flatteningModes = None

    def computeFlatCmd(self, n_modes):
        """
        Compute the command to apply to flatten the input shape.

        Returns
        -------
        flat_cmd : ndarray
            Flat command.
        """
        cmd_amp = -np.dot(self.shape2flat.compressed(), self._recMat)
        _cmd = np.dot(self._cmdMat, cmd_amp)
        #qui impacchettare il flat cmd:
        if isinstance(n_modes, int):
            flat_cmd = self._cmdMat[:,:n_modes] @ _cmd
        elif isinstance(n_modes, list):
            _cmdMat = np.zeros((self._cmdMat.shape[1], len(n_modes)))
            for i,mode in enumerate(n_modes):
                _cmdMat.T[i] = self._cmdMat.T[mode]
            flat_cmd = _cmdMat @ _cmd
        else:
            raise TypeError("n_modes must be either an int or a list of int")
        self.flatCmd = flat_cmd
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
        self.shape2flat = img
        self._rec.loadShape2Flat(img)
        if compute_rec:
            print("Computing recontruction matrix...")
            self._recMat = self._rec.run()

    def reload_intCube(self, tn, zernModes:list=None):
        """
        Reload function for the interaction cube

        Parameters
        ----------
        tn : str
            Tracking number of the new data.
        zernModes : list, optional
            Zernike modes to filter out this cube (if it's not already filtered).
            Default modes are [1,2,3] -> piston/tip/tilt.
        """
        self._intCube = self._loadIntCube(tn, zernModes)
        self._cmdMat  = self._loadCmdMat()
        self._rec     = self._loadReconstructor(self._intCube)
        if self.shape2flat is not None:
            self._recMat = self._rec.run()

    def _loadIntCube(self, zernModes:list=None):
        """
        Interaction cube loader
        
        Parameters
        ----------
        zernModes : list
            Zernike modes to filter out this cube (if it's not already filtered).
            Default modes are [1,2,3] -> piston/tip/tilt.
        
        Return
        ------
        intCube : ndarray
            The interaction cube data array.
        """
        with open(os.path.join(self._path, ifp.flagFile), 'r', encoding='utf-8') as f:
            flag = f.read()
        if ' filtered ' in flag:
            intCube = rd.read_phasemap(os.path.join(self._path, ifp.cubeFile))
        else:
            if zernModes is not None:
                intCube, new_tn = ifp.filterZernikeCube(self._tn, zernModes)
            else:
                intCube, new_tn = ifp.filterZernikeCube(self._tn)
            self.__update_tn(new_tn)
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

    def _loadReconstructor(self, intCube):
        """
        Builds the reconstructor object off the input cube

        Parameters
        ----------
        intCube : ndarray
            Interaction cube data array.

        Returns
        -------
        rec : object
            Reconstructor class.
        """
        rec = crec.ComputeReconstructor(self._intCube)
        return rec

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

    def __update_tn(self,tn):
        """
        Updates the tn and cube path if the tn is to change

        Parameters
        ----------
        tn : str
            New tracking number.
        """
        self._tn = tn
        self._path = os.path.join(ifp.intMatFold, self._tn)
