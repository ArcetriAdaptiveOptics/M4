import numpy as np
from m4.devices import deformable_mirror as dm
from utilf import iff_module as iff
from analyzers import compute_reconstructor as crec
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

    def __init__(self,tn, dm, interf):
        """The Constructor"""
        self._tn = tn
        self._intMat = _loadIntMat()
        #self.dm = dm.AdOpticaDM()?? no, non serve importare il DM
        self._cmdMat = _loadCmdMat()
        self._recMat  = None
        self._frameCenter = None
        self._shape2flat = None
        self._flatOffset = None
        self._cavityOffset = None
        self._synthFlat = None
        self._flatResidue = None
        self._flatResult = None
        self._flatteningModes = None
        pass

    def _loadIntMat(self):
        pass

    def _loadCmdMat(self):
        pass

    def _loadFrameCenter(self):
        pass

    def _computeRec(self, nmodes, zernike2remove):
        pass

    def _registerShape(self, shape, ):
        dp = iffp.findFrameOffset(tn,xxx)
        #cannot work. we should create a dedicated function, not necessarily linked to IFF or flattening
        pass

    def flatten(nmodes):
        pass
