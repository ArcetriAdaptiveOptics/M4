"""
Authors
  - C. Selmi: written in 2022
  - L. Oggioni: added some functionalities
"""

import os as _os
_join = _os.path.join
import logging
import numpy as np
from m4.configuration.ott_parameters import M4Parameters
from m4.configuration import folders as conf
from opticalib.analyzer import modeRebinner
from opticalib.ground.osutils import (
    newtn as _ts, 
    load_fits as _lf,
    save_fits as _sf
)



class FakeM4DM:
    """
    HOW TO USE IT::

        from m4.ott_sim.fake_deformable_mirror import FakeM4DM
        dm = FakeM4DM()
        dm.getActsCommand()
    """

    def __init__(self):
        """The constructor"""
        self.mirrorModes = None
        self.nActs = self.getNActs()
        self.cmdHistory = None
        self._actPos = np.zeros(M4Parameters.N_ACT_SEG)
        self._logger = logging.getLogger("FakeM4DM")
        self._conf = _join(conf.MIRROR_FOLDER, conf.mirror_conf)
        self.m4pupil = _lf(
            _join(self._conf, "m4_mech_pupil-bin2.fits")
        )
        self.m4ima = self.m4pupil * 0.0
        self.CapsensGain = np.load(_join(self._conf, "CapsensGain.npy"))
        self._ifmat = _lf(
            _join(self._conf, "if_sect4_rot-bin2.fits")
        )
        self.ifmat = lambda ifMat: self.getNActs() ** (-0.5) * (
            self._ifmat / np.abs(self._ifmat).max()
        )
        self.ifidx = _lf(
            _join(self._conf, "if_idx4_rot-bin2.fits")
        )
        self.vmat = _lf(
            _join(self._conf, "ff_v_matrix.fits")
        )

    def get_shape(self):
        """
        Returns
        -------
        shape: numpy array
            shape of the mirror
        """
        return self._getActsCommand()

    def set_shape(self, command, differential: bool = False):
        """
        Parameters
        ----------
        command: numpy array [NActs]
            command for a segment
        differential: boolean
            if differential is True the command is added to the previous one
            else the command is absolute
        """
        self._setActsCommand(command, differential)

    def uploadCmdHistory(self, cmdHist):
        """
        Parameters
        ----------
        cmdHist: numpy array [NActs, NFrames]
            command history
        """
        self.cmdHistory = cmdHist

    def runCmdHistory(self, interf=None, rebin: int = 1):
        """
        Run command history

        Parameters
        ----------
        interf: object
            Interferometer to acquire and store measurements.
        """
        baseDataPath = conf.OPD_IMAGES_ROOT_FOLDER
        if self.cmdHistory is None:
            raise Exception("No Command History to run!")
        else:
            tn = _ts.now()
            print(f"{tn} - {self.cmdHistory.shape[-1]} images to go.")
            datafold = _join(baseDataPath, tn)
            if not _os.path.exists(datafold):
                _os.mkdir(datafold)
            for i, cmd in enumerate(self.cmdHistory.T):
                print(f"{i+1}/{self.cmdHistory.shape[-1]}", end="\r", flush=True)
                self.set_shape(cmd)
                if interf is not None:
                    img = interf.acquire_phasemap()
                    img = modeRebinner(img, rebin)
                    path = _join(datafold, f"image_{i:05d}.fits")
                    _sf(path, img)
        self.set_shape(np.zeros(self.nActs))
        return tn

    def _setActsCommand(self, command, rel=False):
        """
        Paramenters
        -----------
        command: numpy array [NActs]
            command for a segment

        Other Parameters
        ----------------
        rel: boolean
            if rel is True relative command is used
            else absolute command
        """
        self._actPos = command
        image = self._mirrorCommand(self._actPos, rel)
        return

    def _getActsCommand(self):
        """
        Returns
        -------
        actsPosition: numpy array [Nacts]
            vector containing the segment actuators position
        """
        return self._actPos

    def getNActs(self):
        """
        Returns
        -------
        n_acts: int
            number of dm actuators
        """
        return M4Parameters.N_ACT_SEG

    def setIncreaseToSegActs(self, inc, rel=True):
        """comando di posizione a tutti i pistoni (old act_incr)

        Parameters
        ----------
        inc: float
            increase in meters
        rel: boolean
            relative increase (True) or absolute position (False)
        """
        inc = inc
        comm = np.ones(self.getNActs()) * inc
        self.set_shape(comm, rel)
        return

    def setZerosToSegActs(self):
        """resets piston position (old act_zero)"""
        comm = np.zeros(self.getNActs())
        self.set_shape(comm)
        return

    def setRandomCommandToSegActs(self, ampiezza, rel=False):
        """generates a random distribution of pistons between
        0 and 'amplitude' (old act_random)

        Parameters
        ----------

        ampiezza: float [m]
            maximum amplitude in meter generated by
            the random distribution of pistons

        rel: boolean
            relative increase (True) or absolute position (False)
        """
        ampiezza = ampiezza
        comm = np.random.rand(self.getNActs()) * ampiezza
        self.set_shape(comm, rel)
        return

    def _generateAndSaveCapsensGain(self, file_name):
        """Function to generate and save a random gain
        Parameters
        ----------
        file_name: string
            cap sens gain file name (.npy)

        Returns
        ------
        file_path: string
            cap sens complete file path
        """
        gain = np.ones(self.getNActs()) + np.random.rand(self.getNActs()) * 0.1 - 0.05
        file_path = _join(self._conf, file_name)
        np.save(file_path, gain)
        return file_path

    def _applyCapsensGain(self, command):
        """
        Parameters
        ----------
        command: numpy array [892]
            vector for multiply capsens gain (for real effect)

        Returns
        ------
        comm: numpy array
            vector multiplied
        """
        command = command * self.CapsensGain
        return command

    def _mirrorCommand(self, comm, rel=True):
        """
        Parameters
        ----------
        comm: numpy array [NActs]
            command for a segment

        Other Parameters
        ----------------
        rel: boolean
            if rel is True relative command is used
            else absolute command

        Return
        ------
        self.m4ima: numpy array
            m4 total image
        """
        comm = self._applyCapsensGain(comm)
        if rel == True:
            self.m4ima.flat[self.ifidx] = self.m4ima.flat[self.ifidx] + np.matmul(
                np.transpose(self.ifmat(self._ifmat)), comm
            )
        elif rel == False:
            self.m4ima.flat[self.ifidx] = np.matmul(
                np.transpose(self.ifmat(self._ifmat)), comm
            )
        return self.m4ima
