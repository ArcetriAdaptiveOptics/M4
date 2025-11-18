"""
Authors
  - L. Oggioni: written in 2022
"""

import os
from re import match
import numpy as np
from m4.configuration import folders as fold_name
from astropy.io import fits
from opticalib.ground import osutils as osu


class ActuatorPositionGenerator:
    """
    HOW TO USE IT::

        from m4.ott_sim.actuator_positions_generator import ActuatorPositionGenerator
        conf = '---/MyConf.yaml'
        ott, interf, dm = start.create_ott(conf)
        ap = ActuatorPositionGenerator(ott, interf, dm)
    """

    def __init__(self, ott, interf, dm):
        """The constructor"""
        self._ott = ott
        self._interf = interf
        self._dm = dm

    def define_actuatorMask_interf_pointOfView(self):
        """
        applica le IF zonali (un attuatore alla volta)
        e crea una maschera attorno ad ognuno applicando una treshold

        poi salva in una matrice complessiva 512*512 x 892 sia le immagini
        che le maschere

        CONTROLLARE CHE FUNZIONI SUGLI ATTUATORI AL BORDO!!!
        """
        self.move2petalo(num=6, RM=0)
        dir0 = os.path.join(fold_name.SIMDATACALIB_ROOT_FOLDER, "PositionActuators", osu.newtn())
        os.makedirs(dir0, exist_ok=True)
        self.dir = dir0

        act_pos_matrix = np.zeros([512 * 512, self._dm.getNActs()])
        act_zonal = np.zeros([512 * 512, self._dm.getNActs()])
        amp = (10 * 1e-9) * 1e6

        for x in range(0, self._dm.getNActs()):
            print("step " + str(x + 1) + "/892", end='\r', flush=True)
            comm = np.zeros(self._dm.getNActs())
            comm[x] = amp
            self._dm.set_shape(comm, differential=False)
            ima = self._interf.acquire_map()
            act_zonal[:, x] = ima.data.flatten()
            act_pos_matrix[:, x] = np.ma.masked_where(ima < 5e-9, ima).mask.flatten()

        osu.save_fits(os.path.join(dir0, "PositionActuator_mask.fits"), act_pos_matrix)
        osu.save_fits(os.path.join(dir0, "PositionActuator_data.fits"), act_zonal)
        return act_pos_matrix, act_zonal

    def check_actposition(self, file: str, N: int):
        if not file in ['data','mask']:
            raise ValueError("file must be 'data' or 'mask'")
        file = os.path.join(self.dir, f"PositionActuator_{file}.fits")
        act_pos_matrix = osu.load_fits(file)
        im = act_pos_matrix[:, N - 1].reshape([512, 512])
        return im

    def move2petalo(self, which: int =1, RM_in: bool = False):
        """
        Move the interferometer to a specific section of the DM

        Parameters
        ----------
        num: integer
            number of the target section
        RM: integer
            Reference mirror in(1) or not(0)
        """
        self._ott.parabolaSlider.setPosition(844)
        if RM_in:
            self._ott.referenceMirrorSlider.setPosition(844)
        if not RM_in:
            self._ott.referenceMirrorSlider.setPosition(0)
        self._ott.angleRotator.setPosition(30 + 60 * (which - 1))
        return
