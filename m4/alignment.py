import numpy as _np
from m4 import folders as _fn
from opticalib import typings as _ot
from opticalib import alignment as _al
from os.path import join as _join, exists as _exists
from opticalib.ground.osutils import newtn as _newtn
from tabulate import tabulate as _tbt
from m4.configuration import ott_status as _otts



class OttAligner(_al.Alignment):
    """
    OttAligner is a class that extends the Alignment class from opticalib.
    It is used to perform alignment operations on optical systems.
    """

    def __init__(self, ott: list[_ot.GenericDevice]|_ot.GenericDevice, interf: _ot.InterferometerDevice) -> None:
        """The Initializer"""
        super().__init__(ott, interf)
        self._parabolatn = _al._sc.fitting_surface.split('/')[-2]
        self._txt = _al._logger.txtLogger(_al._os.path.join(_fn.ALIGNMENT_ROOT_FOLDER, 'AlignmentLog.txt'))


    def correct_alignment(
        self,
        modes2correct: _ot.ArrayLike,
        zern2correct: _ot.ArrayLike,
        applycmd: bool = False,
        n_frames: int = 15,
    ) -> str | _ot.ArrayLike:
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
        coeffs_i = self._zern_routine(self._acquire[0](nframes=1))
        f_cmd = super().correct_alignment(modes2correct=modes2correct, zern2correct=zern2correct, apply=False, n_frames=n_frames)
        if applycmd:
            ntn = _newtn()
            print("Applying correction command...")
            self._apply_command(f_cmd)
            coeffs_f = self._zern_routine(self._acquire[0](nframes=n_frames))
            self._write_correction_log(modes2correct, zern2correct, ntn, coeffs_i, coeffs_f)
            dirr = _join(_fn.ALIGN_RESULTS_ROOT_FOLDER, ntn)
            if not _exists(dirr):
                _al._os.makedirs(dirr)
            _otts.save_positions(dirr, self.mdev)
            self.read_positions()
            return
        else:
            self._write_correction_log([-1], zern2correct, coeffs_i, coeffs_i)
        return f_cmd


    def _write_correction_log(self, m2c: list[int], z2c: list[int], newtn: str, ci: list[float], cf: list[float]) -> None:
        """
        Writes the log of the allignment correction applied to the OTT devices.

        Parameters
        ----------
        initpos : list
            List of the starting positions of the devices, as _Command classes.
        """
        cis, cfs, m2cs, z2cs = _ot.array_str_formatter([ci,cf, m2c, z2c])
        logdict = {
            '': ['Calibration', 'Correction', 'Parabola'], # Rows Names
            'TN':[self._calibtn, newtn, self._parabolatn], # Tracking Numbers
            'Zernike Coeffs': [cis, cfs, f'DoF: {m2cs} | Z2C: {z2cs}'],
        }
        string = _tbt(logdict, headers='keys', tablefmt='plain')
        self._txt.log(string)
        self._txt.log("*"*20)
        return
