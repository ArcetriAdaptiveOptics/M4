"""
Author(s):
----------
- Pietro Ferraiuolo : written in 2024

Description
-----------
This module provides the `Alignment` class and related functions for performing
alignment procedures, including calibration and correction.

How to Use it
-------------
This module contains the class `Alignment`, which manages, alone, both the calibration
and the correction of the alignment of the system. The class is initialized with the
mechanical and acquisition devices used for alignment. These devices, which in
the case of the M4 project are the OTT and the interferometer, are passed as arguments
and configured through the `_m4ac.py` configuration file, (for more information,
check the documentation init).

Given an IPython shell with the tower initialized (`pyott -i`):

>>> from m4.utils.alignment import Alignment
>>> align = Alignment(ott, interf)
>>> # At this point the alignment is ready to be calibrated, given the command amplitude
>>> amps = [0,7, 10, 10, 6, 6, 4, 4] # example, verosimilar, amplitudes
>>> align.calibrate_alignment(amps)
>>> [...]
>>> "Ready for Alignment..."

At this point, the calibration is complete and and `InteractionMatrix.FITS`
was created, saved and stored in the Alignment class. It is ready to compute
and apply corrections. 

>>> modes2correct = [3,4] # Reference Mirror DoF
>>> zern2correct = [0,1] # tip $ tilt
>>> align.correct_alignment(modes2correct, zern2correct, apply=True)

If we already have an InteractionMatrix.FITS file, we can load it and apply 
corrections based off the loaded calibration. All to do is to pass a tracking
number to the `correct_alignment` method:

>>> tn_intmat = '20241122_160000' # example, tracking number
>>> align.correct_alignment(modes2correct, zern2correct, tn=tn_intmat, apply=True)

And the alignment is done.


Notes
-----
Note that the calibration process can be done uploading to the class
a calibrated parabola, so that a different algorithm for the Zernike fitting is 
performed. This can be done through the `reload_calibrated_parabola` method.

>>> tn_parabola = '20241122_160000' # example, tracking number
>>> align.reload_calibrated_parabola(tn_parabola) # load the calibrated parabola
"""

import os
import numpy as np
from m4.utils import _m4ac as mac  # to change
from m4.utils.osutils import findTracknum
from m4.ground import zernike as zern, timestamp as ts, geo
from m4.ground.read_data import readFits_data, readFits_maskedImage, saveFits_data
from m4.ground.logger_set_up import *

_tt = ts.Timestamp()


class Alignment:
    """
    Class for the alignment procedure: calibration and correction.

    This class provides methods to perform alignment procedures using mechanical
    and acquisition devices. It handles the initialization of devices, reading
    calibration data, and executing alignment commands.

    Attributes
    ----------
    mdev : object or list
        The mechanical devices used for alignment. Can be either a single object
        which calls more devices or a list of single devices.
    ccd : object
        The acquisition devices used for alignment.
    cmdMat : numpy.ndarray
        The command matrix read from a FITS file, used for alignment commands.
    intMat : numpy.ndarray or None
        The interaction matrix, initialized as None.
    recMat : numpy.ndarray or None
        The reconstruction matrix, initialized as None.

    Methods
    -------
    correct_alignment(modes2correct, zern2correct, tn=None, apply=False, n_frames=15)
        Corrects the alignment of the system based on Zernike coefficients.
    calibrate_alignment(cmdAmp, n_frames=15, template=None, n_repetitions=1, save=True)
        Calibrates the alignment of the system using the provided command amplitude and template.
    read_positions(show=True)
        Reads the current positions of the devices.
    reload_calibrated_parabola(tn)
        Reloads the calibrated parabola from the given tracking number.

    """

    def __init__(self, mechanical_devices, acquisition_devices):
        """
        Initializes the Alignment class with mechanical and acquisition devices.

        Parameters
        ----------
        mechanical_devices : object or list
            The mechanical devices used for alignment. Can be either
            a single object which calls more devices or a list of
            single devices.
        acquisition_devices : object
            The acquisition devices used for alignment.
        """
        self.mdev = mechanical_devices
        self.ccd = acquisition_devices
        self.cmdMat = readFits_data(mac.commandMatrix)
        self.intMat = None
        self.recMat = None
        self._cmdAmp = None
        self._parabola = (
            readFits_data(mac.calibrated_parabola)
            if not mac.calibrated_parabola == ""
            else None
        )
        self._moveFnc = self.__get_callables(self.mdev, mac.devices_move_calls)
        self._readFnc = self.__get_callables(self.mdev, mac.devices_read_calls)
        self._ccdFncs = self.__get_callables(self.ccd, mac.ccd_acquisition)
        self._devName = self.__get_dev_names(mac.names, ndev=len(self._moveFnc))
        self._dof = [
            np.array(dof) if not isinstance(dof, np.ndarray) else dof for dof in mac.dof
        ]
        self._dofTot = (
            mac.cmdDof
            if isinstance(mac.cmdDof, list)
            else [mac.cmdDof] * len(self._moveFnc)
        )
        self._idx = mac.slices
        self._zvec2fit = np.arange(1, 11)
        self._zvec2use = mac.zernike_to_use
        self._template = mac.push_pull_template
        self._readPath = mac.base_read_data_path
        self._writePath = mac.base_write_data_path
        self._txt = txtLogger(mac.log_path.strip('.log')+'Record.txt')
        set_up_logger(mac.log_path, mac.logging_level)

    def correct_alignment(
        self,
        modes2correct,
        zern2correct,
        tn: str = None,
        apply: bool = False,
        n_frames: int = 15,
    ):
        """
        Corrects the alignment of the system based on Zernike coefficients.

        Parameters
        ----------
        modes2correct : array-like
            Indices of the modes to correct.
        zern2correct : array-like
            Indices of the Zernike coefficients to correct.
        tn : str, optional
            Tracking number of the intMat.fits to be used
        apply : bool, optional
            If True, the correction command will be applied to the system.
            If False (default), the correction command will be returned.
        n_frames : int, optional
            Number of frames acquired and averaged the alignment correction. Default is 15.

        Returns
        -------
        numpy.ndarray or str
            If `apply` is False, returns the correction command as a numpy array.
            If `apply` is True, applies the correction command and returns a string
            indicating that the alignment has been corrected along with the current
            positions.

        Notes
        -----
        This method acquires an image, calculates the Zernike coefficients, reads the
        interaction matrix from a FITS file, reduces the interaction matrix and command
        matrix based on the specified modes and Zernike coefficients, creates a
        reconstruction matrix, calculates the reduced command, and either applies the
        correction command or returns it.
        """
        log(f"{self.correct_alignment.__qualname__}")
        image = self._acquire(n_frames)
        initpos = self.read_positions(show=False)
        zernike_coeff = self._zern_routine(image)
        if tn is not None:
            intMat = readFits_data(self._readPath + f"/{tn}/InteractionMatrix.fits")
            self.intMat = intMat
        else:
            try:
                if self.intMat is not None:
                    intMat = self.intMat
                else:
                    raise AttributeError()
            except AttributeError:
                raise AttributeError(
                    "No internal matrix found. Please calibrate the alignment first."
                )
        reduced_intMat = intMat[
            np.ix_(zern2correct, modes2correct)
        ]
        reduced_cmdMat = self.cmdMat[:, modes2correct]
        recMat = self._create_rec_mat(reduced_intMat)
        reduced_cmd = np.dot(recMat, zernike_coeff[zern2correct])
        f_cmd = -np.dot(reduced_cmdMat, reduced_cmd)
        print(f"Resulting Command: {f_cmd}")
        self._write_correction_log(tn, initpos)
        self._txt.log(
            f"DoF & Zern2Corr:          {modes2correct} {zern2correct}\n" + "-" * 30
        )
        if apply:
            print("Applying correction command...")
            self._apply_command(f_cmd)
            print("Alignment Corrected\n")
            self.read_positions()
            return
        return f_cmd

    def calibrate_alignment(
        self,
        cmdAmp,
        n_frames: int = 15,
        template: list = None,
        n_repetitions: int = 1,
        save: bool = True,
    ):
        """
        Calibrate the alignment of the system using the provided command amplitude and template.

        Parameters
        ----------
        cmdAmp : float
            The command amplitude used for calibration.
        n_frames : int, optional
            The number of frames acquired and averaged for each image. Default is 15.
        template : list, optional
            A list representing the template for calibration. If not provided, the default template will be used.
        n_repetitions : int, optional
            The number of repetitions for the calibration process. Default is 1.
        save : bool, optional
            If True, the resulting internal matrix will be saved to a FITS file. Default is True.

        Returns
        -------
        str
            A message indicating that the system is ready for alignment.

        Notes
        -----
        This method performs the following steps:
        1. Sets the command amplitude.
        2. Uses the provided template or the default template if none is provided.
        3. Produces a list of images based on the template and number of repetitions.
        4. Executes a Zernike routine on the image list to generate an internal matrix.
        5. Optionally saves the internal matrix to a FITS file.
        """
        log(f"{self.calibrate_alignment.__qualname__}")
        self._cmdAmp = cmdAmp
        template = template if template is not None else self._template
        imglist = self._images_production(template, n_frames, n_repetitions)
        intMat = self._zern_routine(imglist)
        self.intMat = intMat
        if save:
            tn = _tt.now()
            filename = os.path.join(self._writePath, tn, "InteractionMatrix.fits")
            os.mkdir(filename.strip("InteractionMatrix.fits"))
            saveFits_data(filename, self.intMat, overwrite=True)
            log(f"{saveFits_data.__qualname__}")
        return "Ready for Alignment..."

    def read_positions(self, show: bool = True):
        """
        Reads the current positions of the devices.

        Returns
        -------
        pos : list
            The list of current positions of the devices.
        """
        log(f"{self.read_positions.__qualname__}")
        logMsg = ""
        pos = []
        logMsg += "Current Positions\n"
        for fnc, dev_name in zip(self._readFnc, self._devName):
            temp = fnc()
            pos.append(_Command(temp))
            logMsg += f"{dev_name}" + " " * (16 - len(dev_name)) + f" : {temp}\n"
        logMsg += "-" * 30
        if show:
            print(logMsg)  #!!! debug only
        return pos

    def reload_calibrated_parabola(self, tn):
        """
        Reloads the parabola from the given file path.

        Parameters
        ----------
        tn : str
            The tracking number of the parabola file.

        Returns
        -------
        str
            A message indicating the successful loading of the file.
        """
        filepath = findTracknum(tn, complete_path=True)
        self._parabola = readFits_maskedImage(filepath)
        return f"Correctly loaded '{filepath}'"

    def _images_production(self, template, n_frames, n_repetitions):
        """
        Produces images based on the provided template and number of repetitions.

        Parameters
        ----------
        template : list
            The template used for image production.
        n_repetitions : int
            The number of repetitions for image production.

        Returns
        -------
        n_results : list
            The list of produced images.
        """
        results = []
        n_results = []
        for i in range(n_repetitions):
            logMsg = ""
            logMsg += f"Repetition n.{i}\n"
            # logging.info(logMsg)
            print(logMsg)  ## debug only
            for k in range(self.cmdMat.shape[1]):
                logMsg2 = ""
                logMsg2 += f"Matrix Column {k+1} : {self.cmdMat.T[k]}"
                print(f"Matrix Column {k+1} : {self.cmdMat.T[k]}\n")  ## debug only
                # logging.info(logMsg2)
                imglist = self._img_acquisition(k, template, n_frames)
                image = self._push_pull_redux(imglist, template)
                image = image / self._cmdAmp[k]
                results.append(image)
            if n_repetitions != 1:
                n_results.append(results)
            else:
                n_results = results
        return n_results

    def _zern_routine(self, imglist):
        """
        Creates the interaction matrix from the provided image list.

        Parameters
        ----------
        imglist : list
            The list of images used to create the interaction matrix.

        Returns
        -------
        intMat : ndarray
            The interaction matrix created from the images.
        """
        log(f"{self._zern_routine.__qualname__}")
        coefflist = []
        if not isinstance(imglist, list):
            imglist = [imglist]
        for img in imglist:
            if self._parabola is None:
                coeff, _ = zern.zernikeFit(img, self._zvec2fit)
                log(f"{zern.zernikeFit.__qualname__}")
            else:
                print("Removing the PAR")
                img = img - 2 * self._parabola
                cir = geo.qpupil(-1 * self._parabola.mask + 1)
                mm = geo.draw_mask(
                    self._parabola.data * 0, cir[0], cir[1], 1.44 / 0.00076 / 2, out=0
                )
                coeff, _ = zern.zernikeFitAuxmask(img, mm, self._zvec2fit)
                log(f"{zern.zernikeFitAuxmask.__qualname__}")
            coefflist.append(coeff[self._zvec2use])
        if len(coefflist) == 1:
            coefflist = [c for c in coefflist[0]]
            return np.array(coefflist)
        intMat = np.array(coefflist).T
        return intMat

    def _create_rec_mat(self, intMat):
        """
        Creates the reconstruction matrix off the inversion of the interaction
        matrix obtained in the alignment calibration procedure.

        Returns
        -------
        recMat : ndarray
            Reconstruction matrix.
        """
        log(f"{self._create_rec_mat.__qualname__}")
        recMat = np.linalg.pinv(intMat)
        self.recMat = recMat
        return recMat

    def _apply_command(self, fullCmd):
        """
        Applies the full command to the devices.

        Parameters
        ----------
        fullCmd : list or ndarray
            Full command of the interaction matrix which commands all device's
            available motors.

        Returns
        -------
        None
        """
        device_commands = self._extract_cmds_to_apply(fullCmd)
        logMsg = ""  #!!!
        for cmd, fnc, dev in zip(device_commands, self._moveFnc, self._devName):
            if cmd.to_ignore:
                logMsg += f"Skipping null command for {dev}\n"  # debug
            else:
                try:
                    logMsg += f"Commanding {cmd} to {dev}\n"  # debug
                    fnc(cmd.vect)
                    log(f"{fnc.__qualname__} : {cmd.vect}")
                except Exception as e:
                    print(e)
                    # logging.warning(f"Someting went wrong with {dev}: {e}")
        logMsg += "-" * 30  # debug
        # logging.info(logMsg)
        print(logMsg)  #!!! debug only

    def _extract_cmds_to_apply(self, fullCmd):
        """
        Extracts the commands to be applied from the full command.

        Parameters
        ----------
        fullCmd : list or ndarray
            The full command from which individual device commands are extracted.

        Returns
        -------
        device_commands : list
            The list of commands to be applied to each device.
        """
        log(f"{self._extract_cmds_to_apply.__qualname__}")
        commands = []
        for d, dof in enumerate(self._dof):
            dev_cmd = np.zeros(self._dofTot[d])
            dev_idx = fullCmd[self._idx[d]]
            for i, idx in enumerate(dev_idx):
                dev_cmd[dof[i]] = idx
            commands.append(_Command(dev_cmd))
        positions = self.read_positions(show=False)
        device_commands = []
        for pos, cmd in zip(positions, commands):
            res_cmd = pos + cmd
            device_commands.append(res_cmd)
        return device_commands

    def _img_acquisition(self, k, template, n_frames):
        """
        Acquires images based on the given template.

        Parameters
        ----------
        k : int
            The index of the command matrix column.
        template : list
            The template used for image acquisition.

        Returns
        -------
        imglist : list
            The list of acquired images.
        """
        log(f"{self._img_acquisition.__qualname__}")
        log(f"{self._acquire.__qualname__}")
        imglist = [self._acquire(n_frames)]
        for t in template:
            logMsg = ""
            logMsg += f"t = {t}"
            cmd = self.cmdMat.T[k] * self._cmdAmp[k] * t
            logMsg += f" - Full Command : {cmd}"
            # logging.info(logMsg)
            print(logMsg)  #!!! debug only
            self._apply_command(cmd)
            log(f"{self._acquire.__qualname__}")
            imglist.append(self._acquire(n_frames))
        return imglist

    def _acquire(self, n_frames):
        """
        Acquires images from the CCD device.

        Parameters
        ----------
        n_frames : int
            The number of frames to acquire and average.

        Returns
        -------
        img : np.ndarray
            The acquired image, correctly reframed according to the acquisition
            device configuration settings.
        """
        return self._ccdFncs[1](self._ccdFncs[0](n_frames))

    def _push_pull_redux(self, imglist, template):
        """
        Reduces the push-pull images based on the given template.

        Parameters
        ----------
        imglist : list
            The list of images to be reduced.
        template : list
            The template used for image reduction.

        Returns
        -------
        image : np.ndarray
            The reduced image.
        """
        log(f"{self._push_pull_redux.__qualname__}")
        template.insert(0, 1)
        image = np.zeros((imglist[0].shape[0], imglist[0].shape[1]))
        for x in range(1, len(imglist)):
            opd2add = imglist[x] * template[x] + imglist[x - 1] * template[x - 1]
            mask2add = np.ma.mask_or(imglist[x].mask, imglist[x - 1].mask)
            if x == 1:
                master_mask = mask2add
            else:
                master_mask = np.ma.mask_or(master_mask, mask2add)
            image += opd2add
        image = np.ma.masked_array(image, mask=master_mask) / 6
        template.pop(0)
        return image

    def _write_correction_log(self, tn, initpos):
        """
        Writes the log of the allignment correction applied to the OTT devices.

        Parameters
        ----------
        initpos : list
            List of the starting positions of the devices, as _Command classes.
        """
        endpos = self.read_positions(show=False)
        par_i, rm_i, m4_i = initpos
        par_f, rm_f, m4_f = endpos
        self._txt.log(
            "Calib. Trackn & IniPos:  {} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e}".format(
                tn, *par_i.vect, *rm_i.vect, *m4_i.vect
            )
        )
        self._txt.log(
            "Result Trackn & EndPos:  {} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e}".format(
                tn, *par_f.vect, *rm_f.vect, *m4_f.vect
            )
        )
        return

    @staticmethod
    def __get_callables(devices, callables):
        """
        Retrieves callable functions for the given devices and calls from
        a configuration file.

        Parameters
        ----------
        devices : object or list
            The devices for which to retrieve callable functions.
        callables : list of str
            The names of the callable functions to retrieve.

        Returns
        -------
        functions : list of callable
            A list of callable functions for the given devices and calls.
        """
        if not isinstance(devices, list):
            devices = [devices]
        functions = []
        for dev in devices:
            for dev_call in callables:
                obj, *methods = dev_call.split(".")
                call = getattr(dev, obj)
                for method in methods:
                    call = getattr(call, method)
                functions.append(call)
        return functions

    @staticmethod
    def __get_dev_names(names, ndev):
        """
        Retrieves device names for the given names and number of devices.

        Parameters
        ----------
        names : list of str
            The names of the devices.
        ndev : int
            The number of devices.

        Returns
        -------
        names : list of str
            A list of device names.
        """
        dev_names = []
        try:
            for x in names:
                dev_names.append(x)
        except TypeError:
            for x in range(ndev):
                dev_names.append(f"Device {x}")
        return dev_names


class _Command:
    """
    The _Command class represents a command with a vector and a flag indicating
    whether the command should be ignored. It provides methods for initializing
    the command, combining it with other commands, and checking if it is null.

    Attributes:
        vect (np.ndarray): The vector representing the command.
        to_ignore (bool): A flag indicating whether the command should be ignored.

    Methods:
        __init__(vector=None, to_ignore:bool=None):
            Initializes a new instance of the _Command class.
        __repr__():
            Returns a string representation of the _Command instance.
        __str__():
            Returns the string representation of the command vector.
        __add__(other):
            Combines the current command with another _Command instance.
        is_null():
            Determines whether the command is null, i.e., a sequence of zeros.
        _process_command_logic(P, C, S):
            Processes the command logic to determine the to_ignore flag.
    """

    def __init__(self, vector=None, to_ignore: bool = None):
        """
        Initializes a new instance of the _Command class.

        Parameters
        ----------
        vector : list or np.ndarray, optional
            The vector representing the command. If a list is provided, it will
            be converted to a numpy array.
        to_ignore : bool, optional
            A flag indicating whether the command should be ignored.
        """
        self.vect = np.array(vector) if isinstance(vector, list) else vector
        self.to_ignore = to_ignore

    def __repr__(self):
        """
        Returns a string representation of the _Command instance.

        Returns
        -------
        str
            A string representation of the _Command instance.
        """
        if self.to_ignore is not None:
            return f"Command({self.vect}, to_ignore={self.to_ignore})"
        else:
            return f"Command({self.vect},)"

    def __str__(self):
        """
        Returns the string representation of the command vector.

        Returns
        -------
        str
            The string representation of the command vector.
        """
        return self.vect.__str__()

    def __add__(self, other):
        """
        Combines the current command with another _Command instance.

        Parameters
        ----------
        other : _Command
            Another instance of the _Command class.

        Returns
        -------
        _Command
            A new _Command instance with the combined vector and updated
            'to_ignore 'flag.

        Raises
        ------
        NotImplementedError
            If the vectors of the commands are not numpy arrays.
        """
        if not isinstance(other, _Command):
            return NotImplemented
        if not isinstance(self.vect, np.ndarray) and not isinstance(
            other.vect, np.ndarray
        ):
            raise NotImplementedError(
                f"Operation not supported for operands types {type(self.vect)} and {type(other)}"
            )
        combined_vect = self.vect + other.vect
        to_ignore = self._process_command_logic(self, other, combined_vect)
        return _Command(combined_vect, to_ignore)

    @property
    def is_null(self):
        """
        Determines whether the command is null, i.e., a sequence of zeros.

        Returns
        -------
        bool
            True if the command is null, False otherwise.
        """
        return np.all(self.vect == 0)

    def _process_command_logic(self, P, C, S):
        """
        Processes the command logic to determine the to_ignore flag.

        Parameters
        ----------
        P : _Command
            The previous command instance.
        C : _Command
            The current command instance.
        S : np.ndarray
            The sum of the vectors of the previous and current commands.

        Returns
        -------
        bool
            The decision for the to_ignore flag based on the command logic.
        """
        # P = current device position
        # C = received device command
        # S = sum of P and C - command to apply (absolute)
        # _________________________________________________#
        # If S = 0
        if np.all(S == 0):
            # C ≠ 0 and P ≠ 0 → TO_NOT_IGNORE
            if not P.is_null and not C.is_null and np.array_equal(C.vect, -1 * P.vect):
                decision = False
            # C = 0 and P = 0 → TO_IGNORE
            elif C.is_null and P.is_null:
                decision = True
        # If S ≠ 0
        else:
            # P ≠ 0 and C = 0 → TO_IGNORE
            if not P.is_null and C.is_null and np.array_equal(S, P.vect):
                decision = True
            # C ≠ 0 and P ≠ 0 → TO_NOT_IGNORE
            elif not C.is_null and not P.is_null:
                decision = False
            # P = 0 and C ≠ 0 → TO_NOT_IGNORE
            elif P.is_null and not C.is_null and np.array_equal(S, C.vect):
                decision = False
        return decision
