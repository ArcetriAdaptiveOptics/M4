"""
Authors
  - L. Oggioni: written in 2022
"""

import os
import numpy as np

# from m4.configuration import start
# from matplotlib import pyplot as plt
from m4.ground import read_data, tracking_number_folder
from m4.configuration import config_folder_names as fold_name
from astropy.io import fits
from opticalib.ground import newtn as ts

# conf='G:\Il mio Drive\Lavoro_brera\M4\LucaConfig.yaml'
# ott, interf, dm = start.create_ott(conf)


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
        dir0 = fold_name.SIMUL_DATA_CALIB_DM_FOLDER + "\\PositionActuators"
        dir, _ = tracking_number_folder.createFolderToStoreMeasurements(dir0)

        act_pos_matrix = np.zeros([512 * 512, self._dm.getNActs()])
        act_zonal = np.zeros([512 * 512, self._dm.getNActs()])

        for x in range(0, self._dm.getNActs()):
            print("step " + str(x + 1) + "/892")
            ampiezza = (10 * 1e-9) * 1e6
            comm = np.zeros(self._dm.getNActs())
            comm[x] = ampiezza
            rel = False
            self._dm.setActsCommand(comm, rel)
            indet = False
            ima = self._interf.acquire_phasemap(indet)

            act_zonal[:, x] = ima.data.flatten()

            im = np.ma.masked_where(ima < 5e-9, ima)
            act_mask = im.mask.flatten()
            act_pos_matrix[:, x] = act_mask

        hdu = fits.PrimaryHDU(act_pos_matrix)
        hdu.writeto(os.path.join(dir, "PositionActuator_mask.fits"))
        hdu = fits.PrimaryHDU(act_zonal)
        hdu.writeto(os.path.join(dir, "PositionActuator_data.fits"))
        return act_pos_matrix, act_zonal

    def check_actposition(self, directory, cartella, N, arg=1):
        name = cartella
        if arg == 1:
            dir = os.path.join(directory, name, "PositionActuator_data.fits")
        if arg == 0:
            dir = os.path.join(directory, name, "PositionActuator_mask.fits")
        act_pos_matrix = read_data.readFits_data(dir)
        im = act_pos_matrix[:, N - 1].reshape([512, 512])
        return im

    def move2petalo(self, num=1, RM=0):
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
        if RM == 1:
            self._ott.referenceMirrorSlider.setPosition(844)
        if RM == 0:
            self._ott.referenceMirrorSlider.setPosition(0)
        self._ott.angleRotator.setPosition(30 + 60 * (num - 1))
        return
