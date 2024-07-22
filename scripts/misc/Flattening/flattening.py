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
        self._intMAt = _loadIntMat()
        #self.dm = dm.AdOpticaDM()??
        pass

    def _loadIntMat(self):
        xxxreadIFF(_self.tn)

    def _computeRec(self, nmodes, zernike2remove):
        pass


