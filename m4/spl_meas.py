import os as _os
import numpy as _np
from tqdm import trange
from m4 import folders as _fn
from opticalib import typings as _ot
from opticalib.ground import osutils as _osu
from photutils import centroids as _centroids
from opticalib.devices.cameras import AVTCamera as _cam
from opticalib.core.fitsarray import fits_array as _fits_array
from opticalib.ground.logger import getSystemLogger as _getLogger


class SplAcquirer:

    def __init__(self, filter: object, camera: str | _cam | None = None):
        """The Constructor"""
        self._filter = filter
        if camera is None:
            camera = "SplCam0"
        elif isinstance(camera, str):
            camera = _cam(name=camera)
        self._camera = camera
        self._logger = _getLogger()

    def acquire(
        self,
        exptime: float,
        lambda_vector: _ot.ArrayLike | None = None,
        numframes: int = 1,
        mask: _ot.MaskData | None = None,
    ):
        """
        Acquire SPL measurements at different wavelengths, by moving the
        "tunable filter".

        Parameters
        ----------
        exptime: float
            Base exposure time of the camera in milliseconds
        lambda_vector : ArrayLike, optional
            Wavelenghts vector, of wavelengths between 400 and 700 nm. If None,
            a default vector is used:
            - from 400 to 700 with 20 nm step

            By default None.
        nframes : int, optional
            number of frames to average for each wavelength, by default 1
        mask: MaskData
            Mask to apply to the measurements. By default, an ampty mask
            is applied.

        Returns
        -------
        tt: string
            Tracking number of measurements, in the `` folder.
        """
        if lambda_vector is not None:
            lambda_vector = _np.asarray(lambda_vector)
            if _np.min(lambda_vector) < 400 or _np.max(lambda_vector) > 700:
                self._logger.error(
                    f"[{self.__class__.__qualname__}] AcquisitionError: Wavelengths must be between 400 and 700 nm"
                )
                raise ValueError("Wavelengths must be between 400 and 700 nm")
        else:
            lambda_vector = _np.arange(400, 701, 20)

        tn = _osu.newtn()
        _osu.save_fits(
            _os.path.join(_fn.OPD_IMAGES_ROOT_FOLDER, tn, "lambda_vector.fits"),
            lambda_vector,
        )
        self._logger.info(
            f"[{self.__class__.__qualname__}] Starting SPL acquisition with tracking number: {tn}"
        )

        ## find PSF position ##
        self._logger.info(
            f"[{self.__class__.__qualname__}] Preparing reference image acquisition"
        )
        self._filter.move_to(600)
        self._camera.set_exptime(exptime / 2 * 1e3)
        self._logger.info(
            f"[{self.__class__.__qualname__}] Acquiring reference image at 600 nm with exposure time of {self._camera.get_exptime()/1000} [ms]"
        )
        img = self._camera.acquire_frames()
        if mask is None:
            mask = _np.zeros(img.shape)
        reference_image = _np.ma.masked_array(img, mask)
        if _np.max(reference_image) > 4000:
            self._logger.warning(
                f"[{self.__class__.__qualname__}] {_np.max(reference_image)} > 4000 : Saturation detected!"
            )

        # Find barycenter
        cy, cx = self._baricenterCalculator(reference_image)

        # Create the Gain Vector
        expgain = _np.ones(lambda_vector.shape[0]) * 0.5
        expgain[_np.where(lambda_vector < 550)] = 1  # 8
        expgain[_np.where(lambda_vector < 530)] = 2  # 8
        expgain[_np.where(lambda_vector > 650)] = 1  # 3
        expgain[_np.where(lambda_vector > 700)] = 1.5  # 8
        self._logger.info(f"[{self.__class__.__qualname__}] Acquisition of frames")

        for wl, expg in zip(lambda_vector, expgain):
            self._logger.info(
                f"Acquiring image at {wl} [nm] with exposure time of {self._camera.get_exptime()/1000} [ms]"
            )
            self._filter.move_to(wl)
            self._camera.set_exptime(exptime * expg * 1e3)
            image = _fits_array(
                data=_np.mean(self._camera.acquire_frames(numframes), 2),
                mask=mask,
                header={
                    "EXPTIME": self._camera.get_exptime() / 1000,
                    "WAVELEN": wl,
                    "EXPGAIN": expg,
                },
            )
            crop = self._preProcessing(image, cy, cx)
            # crop_rot = scin.rotate(crop, 23,  reshape=False)
            crop.writeto(f"image_{wl}nm.fits", overwrite=True)

        self._filter.move_to(600)
        self._logger.info(
            f"[{self.__class__.__qualname__}] Saved tracking number: {tn}"
        )

        return tn

    def _baricenterCalculator(self, reference_image: _ot.ImageData):
        """
        Finds the centroid of the PSF in the reference image

        Parameters:
        -----------
        reference_image : ImageData
            The reference image from which to calculate the centroid.

        Returns:
        --------
        cy, cx : int, int
            The y and x coordinates of the centroid.
        """
        counts, bin_edges = _np.histogram(reference_image, bins=100)
        bin_edges = bin_edges[1:]
        thr = 5 * bin_edges[_np.where(counts == max(counts))]
        img = reference_image.copy()
        cx, cy = _centroids.centroid_com(img, mask=(img < thr))
        return _np.int(cy), _np.int(cx)

    def _preProcessing(self, image: _ot.ImageData, cy: int, cx: int) -> _ot.FitsData:
        """
        Pre-processes the acquired image by subtracting the background and cropping around the PSF.

        Parameters:
        -----------
        image : ImageData
            The acquired image to be pre-processed.
        cy : int
            The y coordinate of the PSF centroid.
        cx : int
            The x coordinate of the PSF centroid.

        Returns:
        --------
        crop : ImageData
            The pre-processed cropped image.
        """
        # FIXME
        xcrop = 145  # 150
        ycrop = 95  # 100

        tmp = _np.zeros((image.shape[0], image.shape[1]))
        tmp[cy - ycrop : cy + ycrop, cx - xcrop : cx + xcrop] = 1
        id_bkg = _np.where(tmp == 0)
        bkg = _np.ma.mean(image[id_bkg])
        img: _ot.FitsData = image - bkg
        crop: _ot.FitsData = img[cy - ycrop : cy + ycrop, cx - xcrop : cx + xcrop]
        return crop


class SplAnalyzer:
    """
    Class used to analyze the images at different wavelengths acquired with the
    SPL camera and Tunable filter system.

    Parameters
    ----------
    tn: str, optional
        Tracking number of the fringes to use in the analysis comparison. By default None.
        It can be set later with the method `setFringesTn`.
    """

    def __init__(self, tn: str | None = None):
        """The constructor"""
        self._logger = _getLogger()
        self.tn_fringes = self.setFringesTn(tn) if tn is not None else None
        self._Qm = None
        self._QmSmooth = None
        self._matrix = None
        self._matrixSmooth = None

    def setFringesTn(self, tn: str):
        """
        Set the tracking number of the fringes to use in the analysis comparison.

        Parameters
        ----------
        tn: string
            tracking number of the fringes to use
        """
        self.tn_fringes = tn
        self._fringes_fold = _os.path.join(_fn.SPL_FRINGES_ROOT_FOLDER, tn)

    def analyzer(self, tn: str) -> tuple[int, int]:
        """
        Analyze measurement data and compare it with synthetic data.

        Parameters
        ----------
        tn: str
            tracking number of the measurement data

        Returns
        -------
        piston: int
            Piston value
        piston_smooth: int
            Piston value after smoothing data
        """

        ntn = _osu.create_data_folder(basepath=_fn.SPL_ROOT_FOLDER)
        self._logger.info(
            f"[{self.__class__.__qualname__}] Analysis of tn = {tn} started."
        )

        lambda_vector = _osu.load_fits(_os.path.join(tn, "lambda_vector.fits"))

        cube, cube_normalized = self.readMeasurement(tn)
        matrix, matrix_smooth = self.matrix_calc(lambda_vector, cube, cube_normalized)
        piston, piston_smooth = self._templateComparison(
            matrix, matrix_smooth, lambda_vector
        )

        self._savePistonResult(ntn, piston, piston_smooth)
        # FIXME: save matrix in the results folder
        _osu.save_fits(_os.path.join(ntn, "fringe_result.fits"), matrix)

        return piston, piston_smooth

    def matrix_calc(
        self,
        lambda_vector: _ot.ArrayLike,
        cube: _ot.CubeData,
        cube_normalized: _ot.CubeData,
    ) -> tuple[_ot.MatrixLike, _ot.MatrixLike]:
        """
        Calculate the matrix of fringes from the acquired cube of images.

        Parameters
        ----------
        lambda_vector: numpy array
            Vector of wavelengths (between 400/700 nm)
        cube: numpy array
            Cube of images [pixels, pixels, n_frames=lambda]
        cube_normalized: numpy array
            Cube of normalized images [pixels, pixels, n_frames=lambda]

        Returns
        -------
        matrix: numpy array
            Matrix of fringes
        matrix_smooth: numpy array
            Smoothed matrix of fringes
        """
        img = _np.sum(cube_normalized, 2)
        pick = self._newThr(img)
        matrix = _np.zeros(
            (pick[3] - pick[2] + 1, lambda_vector.shape[0])
        )  # 150 + 1 pixel
        matrix_smooth = _np.zeros((pick[3] - pick[2] + 1, lambda_vector.shape[0]))
        crop_frame_cube = None

        for i in range(lambda_vector.shape[0]):
            frame = cube[:, :, i]
            crop_frame = frame[pick[0] : pick[1], pick[2] : pick[3] + 1]

            if crop_frame_cube is None:
                crop_frame_cube = crop_frame
            else:
                crop_frame_cube = _np.dstack((crop_frame_cube, crop_frame))

            y = _np.sum(crop_frame, 0)
            area = _np.sum(y[:])
            y_norm = y / area
            matrix[:, i] = y_norm

            w = self.smooth(y_norm, 4)
            w = w[: pick[3] - pick[2] + 1]
            matrix_smooth[:, i] = w

        matrix[_np.where(matrix == _np.nan)] = 0
        self._matrix = matrix
        self._matrixSmooth = matrix_smooth
        return matrix, matrix_smooth

    def smooth(self, a: _ot.ArrayLike, WSZ: int):
        """'

        Parameters
        ----------
        a: NumPy 1-D array
            containing the data to be smoothed
        WSZ: int
            smoothing window size needs, which must be odd number,
            as in the original MATLAB implementation

        Returns
        -------
        smooth: numpy array
                smoothd data
        """
        out0 = _np.convolve(a, _np.ones(WSZ, dtype=int), "valid") / WSZ
        r = _np.arange(1, WSZ - 1, 2)
        start = _np.cumsum(a[: WSZ - 1])[::2] / r
        stop = (_np.cumsum(a[:-WSZ:-1])[::2] / r)[::-1]
        return _np.concatenate((start, out0, stop))

    def get_matrix(self, tn: str) -> _ot.MatrixLike:
        """
        Loads the matrix of fringes from a specific tracking number.

        Paramenters
        -----------
        tn: string
           Ttracking number from which to read matrix

        Returns
        -------
        matrix: numpy array
            Matrix of fringes
        """
        return _osu.load_fits(
            _os.path.join(_fn.SPL_ROOT_FOLDER, tn, "fringe_result.fits")
        )

    def readMeasurement(self, tn: str) -> tuple[_ot.CubeData, _ot.CubeData]:
        """
        Read images in a specific tracking number and return the cube of images
        and the cube of normalized images.

        Parameters
        ----------
        tn: string
            tracking number of the measurement data

        Returns
        -------
        cube: numpy array
            Cube of images [pixels, pixels, n_frames=lambda]
        cube_normalized: numpy array
            Cube of normalized images [pixels, pixels, n_frames=lambda]
        """
        cube = _osu.loadCubeFromFilelist(tn, fold=_fn.SPL_ROOT_FOLDER, key="image")
        cube_normalized = map(
            lambda x: x / _np.sum(x), cube.transpose(2, 0, 1)
        ).transpose(1, 2, 0)
        return cube, cube_normalized

    def _newThr(self, img: _ot.ImageData) -> list[int]:
        """
        Calculate the peak position of the image

        Parameters
        ----------
        img: numpy array
            2D image

        Returns
        -------
        peak: list
            list of 4 integers defining the crop area
        """
        cx, cy = _centroids.centroid_2dg(img)
        baricenterCoord = [_np.int(round(cy)), _np.int(round(cx))]
        peak = [
            baricenterCoord[0] - 25,
            baricenterCoord[0] + 25,
            baricenterCoord[1] - 75,
            baricenterCoord[1] + 75,
        ]
        return peak

    def _templateComparison(
        self,
        matrix: _ot.MatrixLike,
        matrix_smooth: _ot.MatrixLike,
        lambda_vector: _ot.ArrayLike,
    ) -> tuple[int, int]:
        """
        Compare the matrix obtained from the measurements with
        the one recreated with the synthetic data in tn_fringes.

        Parameters
        ----------
        matrix: MatrixLike
            Measured matrix, [pixels, lambda]
        matrix_smooth: MatrixLike
            Measured smoothed matrix, [pixels, lambda]
        lambda_vector: ArrayLike
            Vector of wavelengths
        Returns
        -------
        piston: int
                piston value
        """
        self._logger.debug(
            f"[{self.__class__.__qualname__}] Template Comparison with data in {self.tn_fringes}"
        )
        delta, lambda_synth = self._readDeltaAndLambdaFromFringesFolder()
        idx = _np.isin(lambda_synth, lambda_vector)

        Qm = matrix - _np.mean(matrix)
        Qm_smooth = matrix_smooth - _np.mean(matrix_smooth)
        self._Qm = Qm
        self._QmSmooth = Qm_smooth

        F = []
        for i in range(1, delta.shape[0]):
            file_name = _os.path.join(self._fringes_fold, "Fringe_%05d.fits" % i)
            fringe = _osu.load_fits(file_name)
            fringe_selected = fringe[:, idx]
            F.append(fringe_selected)
        F = _np.dstack(F)
        Qt = F - _np.mean(F)
        self._Qt = Qt

        R = _np.zeros(delta.shape[0] - 1)
        R_smooth = _np.zeros(delta.shape[0] - 1)
        for i in trange(delta.shape[0] - 1, desc=f"Comparing with synthetic data"):

            R[i] = _np.sum(Qm[:, :] * Qt[:, :, i]) / (
                _np.sum(Qm[:, :] ** 2) ** 0.5 * _np.sum(Qt[:, :, i] ** 2) ** 0.5
            )
            R_smooth[i] = _np.sum(Qm_smooth[:, :] * Qt[:, :, i]) / (
                _np.sum(Qm_smooth[:, :] ** 2) ** 0.5 * _np.sum(Qt[:, :, i] ** 2) ** 0.5
            )

        idp = _np.where(R == max(R))
        idp_smooth = _np.where(R_smooth == max(R_smooth))
        piston = delta[idp]
        piston_smooth = delta[idp_smooth]
        return piston, piston_smooth

    def _readDeltaAndLambdaFromFringesFolder(
        self,
    ) -> tuple[_ot.ArrayLike, _ot.ArrayLike]:
        """
        Reads the delta piston and synthetic wavelength data from the fringes folder.

        Returns
        -------
        delta: ArrayLike
            Delta piston values
        lambda_synth_from_data: ArrayLike
            Synthetic wavelength values
        """
        delta = _osu.load_fits(
            _os.path.join(self._fringes_fold, "Differential_piston.fits")
        )
        lambda_synth_from_data = (
            _osu.load_fits(_os.path.join(self._fringes_fold, "Lambda.fits")) * 1e9
        ).astype(int)
        return delta, lambda_synth_from_data

    def _savePistonResult(self, tn: str, piston: _ot.Number, piston_smooth: _ot.Number):
        """
        Save the piston results to a text file.

        Parameters
        ----------
        tn: string
            tracking number where to save the results
        piston: int
            piston value
        piston_smooth: int
            piston value after smoothing data
        """
        savepath = _os.path.join(tn, "pistons.txt")
        with open(savepath, "w+") as file:
            file.write(f"[{tn.split('/')[-1]} {piston[0]:.4e}, {piston_smooth[0]:.4e}]")
