import os as _os
import numpy as _np
from tqdm import trange
from m4 import folders as _fn
from opticalib import typings as _ot
from opticalib.ground import osutils as _osu
from photutils import centroids as _centroids
from opticalib.devices.cameras import AVTCamera as _cam
from opticalib.ground.logger import SystemLogger as _SL
from opticalib.core.fitsarray import fits_array as _fits_array

def _get_tunable_filter():
    """
    initiate the tunable filter with standard parameters
    """
    from plico_motor import motor # type: ignore
    from opticalib.core.read_config import getInterfConfig
    
    ip = getInterfConfig('PhaseCam6110')['ip']
    return motor(ip, 7200, axis=0)

class SplAcquirer:

    def __init__(self, camera: str | _cam | None = None, tunable_filter: object | None = None):
        """The Constructor"""
        if camera is None:
            camera = "SplCam0"
        elif isinstance(camera, str):
            camera = _cam(name=camera)
        if tunable_filter is None:
            tunable_filter = _get_tunable_filter()
        self._camera = camera
        self._filter = tunable_filter
        self._darkFrame1sec= None
        self._curr_exptime = None
        self._logger = _SL(__class__)

    def set_exptime(self, exptime: float):
        """
        Set the SPL camera exposure time, if the target value is different from the current one

        Parameters:
        -----------
        exptime : float
            the exposure time in [s]

        Returns:
        --------
        """
        if self._curr_exptime != exptime:
            self._camera.set_exptime(exptime*1e6)
            self._curr_exptime = exptime
        else:
            self._logger.warning("The requested exposure time for the camera is equal to the current one. Skipping")
            pass


    def acquireDarkFrame(self, exptime: float, nframes: int = 1):
        """
        Calibrates a dark frame (or better a sky frame) to allow the subtraction
        of the daylight signal.

        Parameters:
        -----------
        exptime : float
            the exposure time in [s]
        nframes : int
            the number of frames to be averaged together

        Returns:
        --------
        the dark frames scaled to 1sec exposure time
        """
        self.set_exptime(exptime)
        self._logger.info(f'Acquiring Dark frame for {exptime} s exposure time, averaging {nframes} frames')
        img = self._camera.acquire_frames(nframes)
        if not hasattr(img, 'mask'):
            img = _np.ma.masked_array(img, mask=_np.zeros_like(img))
        self._darkFrame1sec = img/exptime
        return self._darkFrame1sec.copy()

    def _removeDarkFrame(self, img: _ot.ImageData, exptime: float) -> _ot.ImageData:
        """
        Corrects a frame for the dark frame and (when implemented) for sky and RON

        Parameters:
        -----------
        img : ImageData
            The frame yto which subtract the Dark
        exptime : float
            The exposure time in s

        Returns:
        --------
        the corrected image
        """

        if self._darkFrame1sec is not None:
            imgout = img - self._darkFrame1sec * exptime
            self._logger.info(f'Subtracting the dark frame scaled for {exptime} s exposure time')
        else:
            imgout = img
        return imgout

    def acquire(
        self,
        exptime: float,
        lambda_vector: _ot.ArrayLike | None = None,
        nframes: int = 1,
        mask: _ot.MaskData | None = None,
    ):
        """
        Acquire SPL measurements at different wavelengths, by moving the
        "tunable filter".

        Parameters
        ----------
        exptime: float
            Base exposure time of the camera in seconds
        lambda_vector : ArrayLike, optional
            Wavelenghts vector, of wavelengths between 400 and 700 nm. If None,
            a default vector is used:
            - from 400 to 700 with 20 nm step

            By default None.
        nframes : int, optional
            number of frames to average for each wavelength, by default 1
        mask: MaskData | None, optional
            Mask to apply to the measurements. By default, an ampty mask
            is applied.

        Returns
        -------
        tn: string
            Tracking number of measurements, in the `` folder.
        """
        if lambda_vector is not None:
            lambda_vector = _np.asarray(lambda_vector)
            if _np.min(lambda_vector) < 400 or _np.max(lambda_vector) > 700:
                self._logger.error(
                    f"AcquisitionError: Wavelengths must be between 400 and 700 nm"
                )
                raise ValueError("Wavelengths must be between 400 and 700 nm")
        else:
            lambda_vector = _np.arange(440, 721, 20)

        datapath = _osu.create_data_folder(basepath=_fn.SPL_DATA_ROOT_FOLDER)
        tn = datapath.split("/")[-1]
        print(tn)

        _osu.save_fits(
            _os.path.join(datapath, "lambda_vector.fits"),
            lambda_vector,
        )

        if self._darkFrame1sec is not None:
            _osu.save_fits(_os.path.join(datapath, 'darkFrame.fits'), self._darkFrame1sec)

        self._logger.info(
            f"Starting SPL acquisition with tracking number: {tn}"
        )

        ## find PSF position ##
        self._logger.info(
            f"Preparing reference image acquisition"
        )
        self._filter.move_to(600)
        self.set_exptime(exptime/2)
        self._logger.info(
            f"Acquiring reference image at 600 nm with exposure time of {self._curr_exptime} [s]"

        )
        img = self._camera.acquire_frames()
        if mask is None:
            mask = _np.zeros(img.shape)
        reference_image = _np.ma.masked_array(img, mask)
        if _np.max(reference_image) > 4000:
            self._logger.warning(
                f"{_np.max(reference_image)} > 4000 : Saturation detected!"
            )

        # Create the Gain Vector
        expgain = _np.ones(lambda_vector.shape[0]) * 0.5
        expgain[_np.where(lambda_vector < 550)] = 1  # 8
        expgain[_np.where(lambda_vector < 530)] = 2  # 8
        expgain[_np.where(lambda_vector > 650)] = 1  # 3
        expgain[_np.where(lambda_vector > 700)] = 1.5  # 8
        self._logger.info(f"Acquisition of frames")

        for wl, t_int in zip(lambda_vector, (expgain*exptime)):
            self._logger.info(
                f"Acquiring image at {wl:.1f} [nm] with exposure time of {self._curr_exptime:.3f} [s]"

            )
            print(f'Moving to lambda: {wl}')
            self._filter.move_to(wl)
            self.set_exptime(t_int)
            img = self._camera.acquire_frames(nframes)
            image = _fits_array(
                data=img,
                mask=mask,
                header={
                    "EXPTIME": self._curr_exptime,
                    "WAVELEN": wl,
                    "EXPGAIN": t_int/self._curr_exptime,
                },
            )
            image.writeto(_os.path.join(datapath, f"rawframe_{wl}nm.fits"), overwrite=True)

        self._filter.move_to(600)
        self._logger.info(
            f"Saved tracking number: {tn}"
        )

        self._last_measure_tn = tn

        return tn

    def postProcessRawFrames(self, tn: str = None, remove_median: bool = True, remove_dark: bool = True, crop: bool = True):
        """
        Takes the raw frames acquired with the `acquire` method and
        applies a soft processing controlled by input parameters.
        
        This processing is applied to every `rawframe_[wavelength]nm.fits` file
        in the last measurement tracking number folder, and the result is saved
        as `postprod_[wavelength]nm.fits`.
        
        Parameters:
        -----------
        remove_median : bool
            If True, removes the median value from each frame.
        remove_dark : bool
            If True, removes the dark frame from each frame.
        """
        tn = self._last_measure_tn if tn is None else tn
        if not hasattr(self, '_last_measure_tn'):
            self._logger.error("No measurement tracking number found. Please run `acquire` method first.")
            raise AttributeError("No measurement tracking number found. Please run `acquire` method first.")

        datapath = _os.path.join(_fn.SPL_DATA_ROOT_FOLDER, tn)
        filelist = _osu.getFileList(
            tn, 
            fold=_fn.SPL_DATA_ROOT_FOLDER, 
            key="rawframe"
        )
        rawlist = _osu.loadCubeFromFilelist(filelist)

        for img, filename in zip(rawlist, filelist):

            if remove_dark:
                img = self._removeDarkFrame(img, exptime=img.header['EXPTIME'])

            if remove_median:
                median_value = _np.median(img)
                img = _np.clip(img - median_value, 0, None)
                
            if crop:
                cy, cx = self._baricenterCalculator(img)
                img = img[cy - 100 : cy + 100, cx - 150 : cx + 150]

            new_filename = filename.replace("rawframe_", "postprod_")
            new_filepath = _os.path.join(datapath, new_filename)
            _osu.save_fits(new_filepath, img)
        
        self._logger.info(
            f"Post-processed tn: {tn}. Removed median: {remove_median}, Removed dark: {remove_dark}"
        )


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
        return int(_np.round(cy)), int(_np.round(cx))


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
        self._logger = _SL(__class__)
        self.tn_fringes = self.setFringesTn(tn) if tn is not None else None
        self._Qm = None
        self._QmSmooth = None
        self._matrix = None
        self._matrixSmooth = None
        self._darkFrame1sec = None


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

    def analyzer(self, tn: str, remove_dark: bool = False, remove_median: bool = False, crop: bool = False) -> tuple[int, int]:
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

        ntn = _osu.create_data_folder(basepath=_fn.SPL_RESULTS_ROOT_FOLDER)
        self._logger.info(
            f"Analysis of tn = {tn} started."
        )

        if any([x is True for x in [remove_dark, remove_median, crop]]):
            self.postProcessRawFrames(tn, remove_median, remove_dark, crop)
            raw = False
        else:
            raw = True

        datapath = _os.path.join(_fn.SPL_DATA_ROOT_FOLDER, tn)

        lambda_vector = _osu.load_fits(_os.path.join(datapath, "lambda_vector.fits" ))

        cube, cube_normalized = self.readMeasurement(tn, raw)
        matrix, matrix_smooth = self.matrix_calc(lambda_vector, cube, cube_normalized)
        piston, piston_smooth = self._templateComparison(
            matrix, matrix_smooth, lambda_vector
        )

        self._savePistonResult(ntn, piston, piston_smooth)
        # FIXME: save matrix in the results folder
        _osu.save_fits(_os.path.join(ntn, "fringe_result.fits"), matrix)

        return piston, piston_smooth

    def fringes(self, tn):
        from matplotlib import pyplot as plt

        cube, _ = self.readMeasurement(tn)
        wav = self.get_measurement_lambda(tn)
        mat, _ = self.matrix_calc(wav,cube,cube)
        npix = mat.shape[0]
        extent2use = [0,npix-1,wav[0],wav[-1]]
        plt.imshow(mat.T,extent = extent2use,aspect =0.3)
        plt.xlabel('Pixel')
        plt.ylabel('Wavelength')
        plt.title(tn)
        return mat
    

    def get_measurement_lambda(self,tn):
        wav = _osu.load_fits(
            _os.path.join(
                _fn.SPL_DATA_ROOT_FOLDER, tn, "lambda_vector.fits"
            )
        )      
        return wav

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

    def postProcessRawFrames(self, tn: str = None, remove_median: bool = True, remove_dark: bool = True, crop: bool = True):
        """
        Takes the raw frames acquired with the `acquire` method and
        applies a soft processing controlled by input parameters.
        
        This processing is applied to every `rawframe_[wavelength]nm.fits` file
        in the last measurement tracking number folder, and the result is saved
        as `postprod_[wavelength]nm.fits`.
        
        Parameters:
        -----------
        remove_median : bool
            If True, removes the median value from each frame.
        remove_dark : bool
            If True, removes the dark frame from each frame.
        """
        datapath = _os.path.join(_fn.SPL_DATA_ROOT_FOLDER, tn)
        filelist = _osu.getFileList(
            fold=datapath,
            key="rawframe"
        )
        rawlist = [_osu.load_fits(x) for x in filelist]

        for img, filename in zip(rawlist, filelist):

            if remove_dark:
                dark = _osu.load_fits(_os.path.join(datapath, 'darkFrame.fits'))
                img = img - dark*img.header['EXPTIME']
                img.header['RDARK'] = True

            if remove_median:
                median_value = _np.ma.median(img)
                #img = (img - median_value).clip(0, None)
                img = img-median_value
                img[img < 0] = 0
                img.header['RMEDIAN'] = True

            if crop:
                cy, cx = self._baricenterCalculator(img)
                img = img[cy - 100 : cy + 100, cx - 150 : cx + 150]
                img.header['CROPPED'] = True

            new_filename = filename.replace("rawframe_", "postprod_")
            new_filepath = _os.path.join(datapath, new_filename)
            img.writeto(new_filepath, overwrite=True)
#            _osu.save_fits(filepath=new_filepath, data=img, overwrite=True, header=header)

        self._logger.info(
            f"Post-processed tn: {tn}. Removed median: {remove_median}, Removed dark: {remove_dark}"
        )


    def readMeasurement(self, tn: str, raw: bool = True) -> tuple[_ot.CubeData, _ot.CubeData]:
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
        if raw:
            key = 'rawframe'
        else:
            key = 'postprod'
        cube = _osu.loadCubeFromFilelist(tn, fold='SPL',  key=key)
        cube = _np.transpose(cube,[1,0,2])  #modRB20250117 to manage now image orientation
        cube_normalized = _np.array(list(map(
            lambda x: x / _np.sum(x), cube.transpose(2, 0, 1)
        ))).transpose(1, 2, 0)
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

        baricenterCoord = [_np.int_(round(cy)), _np.int_(round(cx))]
        peak = [
            baricenterCoord[0] - 25,
            baricenterCoord[0] + 25,
            baricenterCoord[1] - 75,
            baricenterCoord[1] + 75,
        ]
        return peak


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
        return int(_np.round(cy)), int(_np.round(cx))


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
            f"Template Comparison with data in {self.tn_fringes}"
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
        #lambda_synth_from_data = (_osu.load_fits(_os.path.join(self._fringes_fold, "Lambda.fits")) * 1e9  ).astype(int)  #modRB20250117 since astype(int) yields errors
        lambda_synth_from_data = ( _osu.load_fits(_os.path.join(self._fringes_fold, "Lambda.fits")) * 1e9 )
        lambda_synth_from_data = (_np.round(lambda_synth_from_data/5)*5).astype(int)  #end of modification

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
