"""
Author(s) 
---------
    - Pietro Ferraiuolo : written in 2024

Description
-----------
"""
import os
import numpy as np
from m4.devices import deformable_mirror as dm
from m4.dmutils import iff_processing as ifp
from m4.ground import read_data as rd
from m4.utils import osutils as osu
from m4.analyzers import compute_reconstructor as crec
#  ---- leggere flattening.py come possibile riferimento -----

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
        self._recMat        = self._rec.run()
        self._frameCenter   = None
        self._shape2flat    = None # Immagine acquisita o fornita, come? 
        self._flatOffset    = None
        self._cavityOffset  = None
        self._synthFlat     = None
        self._flatResidue   = None
        self._flatResult    = None
        self._flatteningModes = None

    def flatten(self, nmodes):
        pass

    def computeFlatCmd(self):
        # Logica
        cmd_amp = -np.dot(self._recMat.flatten(), self._shape2flat.compressed())
        flat_cmd = np.dot(self._cmdMat, cmd_amp)
        self._flatResult = flat_cmd
        return flat_cmd

    def _loadIntMat(self):
        cube = os.path.join(self._path, ifp.cubeFile)
        with open(os.path.join(self._path, ifp.flagFile), 'r', encoding='utf-8') as f:
            flag = f.read()
        if 'filtered' in flag:
            intmat = rd.read_phasemap(cube)
        else:
            intmat = ifp.filterZernikeCube(self._tn)
        return intmat

    def _loadCmdMat(self):
        cmdMat = rd.readFits_data(os.path.join(self._path, ifp.cmdMatFile))
        return cmdMat

    def _loadFrameCenter(self):
        frame_center = rd.readFits_data('data')
        return frame_center

    def _registerShape(self, shape):
        xxx=None
        dp = ifp.findFrameOffset(self._tn,xxx)
        #cannot work. we should create a dedicated function, not necessarily linked to IFF or flattening
        pass

    def load_image2shape(self, img):
        self._shape2flat = img
        return
