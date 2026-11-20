"""
SPLChain: one tunable filter shared by several single-PSF cameras.

The M4 Sensor for Phase Lag uses ONE tunable filter (Thorlabs Kurios, driven
through plico) shared by SIX GigE cameras, each seeing ONE PSF. The wavelength
is a shared, serial resource, while everything between two filter moves is
independent per camera. `SPLChain` therefore composes one single-camera
``opticalib.procedures.phasing.SPL`` per camera, and overlaps exposure updates
with filter moves, acquires all cameras simultaneously and saves frames in the
background. The filter is only ever driven from the calling thread, and never
while a camera is exposing.

Timing model
------------
For each wavelength: exposure updates (threads) run concurrently with the
filter move (calling thread), then a barrier, then all cameras acquire in
parallel, and the frames are saved by a writer pool while the loop advances to
the next wavelength. For N wavelengths the wall time is about
``N * (T_move + max T_exp)`` instead of ``N * (T_move + 6 * T_exp)`` for a
camera-after-camera sweep.

Data layout
-----------
Under the SPL data root folder::

    SPL/<tn>/lambda_vector.fits
    SPL/<tn>/cam{k}/lambda_vector.fits
    SPL/<tn>/cam{k}/darkFrame.fits
    SPL/<tn>/cam{k}/rawframe_<wl>nm.fits

Frame headers carry ``EXPTIME``, ``WAVELEN``, ``EXPGAIN`` (the gain actually
applied) and ``CAMNAME``. Each ``cam{k}`` folder is a valid single-camera SPL
dataset, analysable on its own with ``SPL``.

Configuration
-------------
Defaults come from the ``PHASING`` section of the configuration file::

    PHASING:
      cameras:                 # one per PSF, in "cameras:<Name>" form
        - "cameras:<Name1>"
        - "cameras:<Name2>"
        - "cameras:<Name3>"
        - "cameras:<Name4>"
        - "cameras:<Name5>"
        - "cameras:<Name6>"
      psfs_angles: [0, 0, 0, 0, 0, 0]   # PSF rotation angle, one per camera
      filter: MOTORS:TunableFilter      # tunable filter, "<DEVICES subsection>:<Name>"

where ``DEVICES.MOTORS.TunableFilter`` holds the ``ip`` and ``port`` of the
plico motor server of the filter.

Requires ``opticalib>=2.5.5``.
"""

import os as _os
import numpy as _np
import multiprocessing as _mp
from concurrent import futures as _fut

from opticalib.core import _types as _ot
from opticalib.core.root import folders as _fn
from opticalib.ground import osutils as _osu
from opticalib.ground.logger import SystemLogger as _SL
from opticalib.core.fitsarray import fits_array as _fits_array
from opticalib.core.config import LazySection as _LazySection
from opticalib.procedures import phasing as _phasing
from opticalib.procedures.phasing import SPL, _read_filter_mode

_splconf = _LazySection("PHASING")


class SPLChainError(RuntimeError):
    """
    Error raised when one or more cameras fail during a chain operation.

    Attributes
    ----------
    errors : dict[int, BaseException]
        Camera index -> exception raised by that camera.
    tn : str | None
        Tracking number of the (possibly partial) measurement, when raised by
        `SPLChain.acquire` after its data folder was created; else None.
    """

    def __init__(
        self,
        message: str,
        errors: dict[int, BaseException],
        tn: str | None = None,
    ):
        super().__init__(message)
        self.errors = errors
        self.tn = tn


def _chain_error(what: str, errors: dict[int, BaseException]) -> SPLChainError:
    """Build a `SPLChainError` summarising the per-camera `errors` of `what`."""
    detail = "; ".join(f"cam{k}: {e!r}" for k, e in sorted(errors.items()))
    return SPLChainError(f"{what} failed on camera(s) {sorted(errors)}: {detail}", errors)


def _add_note(exc: BaseException, note: str) -> None:
    """Attach `note` to `exc` (PEP 678)."""
    exc.add_note(note)


class SPLChain:
    """
    Chain of single-camera SPL sharing one tunable filter.

    Each camera is wrapped by its own ``SPL`` instance (see `spls`); the chain
    coordinates the shared filter and runs the per-camera work in parallel.
    See the module documentation for the timing model and the data layout.

    Failure policy: if a camera fails during a step, the other cameras finish
    that step, then the sweep stops, the filter is parked at 600 nm, the cached
    exposure times are invalidated and `SPLChainError` is raised (``errors``
    keyed by camera index, ``tn`` the tracking number of the partial data).

    The chain can be used as a context manager, which closes the cameras on
    exit.

    Parameters
    ----------
    cameras : list[str | CameraDevice] | None
        Cameras, one per PSF. Strings are camera names (optionally in the
        configuration form ``"cameras:<Name>"``) and are opened as
        `GigaVision`. If None, read from ``PHASING.cameras``.
    tunable_filter : object | None
        The tunable filter client. If None it is created once from the
        configuration. The string ``"none"`` disables it.
    tnfringes : str | None
        Tracking number of the simulated fringes, for template comparison.
    angles : list[float] | None
        One PSF rotation angle per camera. If None, read from
        ``PHASING.psfs_angles``, else zeros.

    Raises
    ------
    ValueError
        If the number of angles differs from the number of cameras. Cameras
        opened by the constructor are closed before any exception propagates.

    Examples
    --------
    .. code-block:: python

        from m4.spl import SPLChain

        with SPLChain() as chain:          # cameras/filter from PHASING
            chain.set_tn_fringes("20240101_120000")
            chain.acquire_dark_frame(exptime=1.0, nframes=5)
            tn = chain.acquire(exptime=0.01, filter_mode="narrow")
            pistons = chain.analysis(tn)   # one value per camera

    With ``max_workers > 1`` (the default) `analysis` spawns processes, so a
    script must be guarded by ``if __name__ == "__main__":``.
    """

    def __init__(
        self,
        cameras: list[str | _ot.CameraDevice] | None = None,
        tunable_filter: object | None = None,
        tnfringes: str | None = None,
        angles: list[float] | None = None,
    ):
        """Open the cameras, the filter and the per-camera SPL instances."""
        self._logger = _SL(__class__)

        if cameras is None:
            cameras = list(_splconf["cameras"])
        opened = []  # cameras created here (not the caller's), closed on failure
        try:
            self._init_members(cameras, tunable_filter, tnfringes, angles, opened)
        except BaseException:
            for cam in opened:
                if hasattr(cam, "close"):
                    try:
                        cam.close()
                    except Exception as e:
                        self._logger.error(f"Error closing camera: {e!r}")
            raise

    def _init_members(self, cameras, tunable_filter, tnfringes, angles, opened):
        cams, names = [], []
        for k, cam in enumerate(cameras):
            if isinstance(cam, str):
                name = cam.split(":")[-1]
                if name.lower() == "none":
                    cam = None
                else:
                    from opticalib.devices.cameras import GigaVision

                    cam = GigaVision(name=name)
                    opened.append(cam)
            else:
                name = getattr(cam, "name", None) or f"cam{k}"
            cams.append(cam)
            names.append(str(name))
        self._names = names

        if angles is None:
            try:
                angles = list(_splconf["psfs_angles"])
            except KeyError:
                angles = [0.0] * len(cams)
        if len(angles) != len(cams):
            raise ValueError(
                f"Got {len(angles)} angles for {len(cams)} cameras: they must match"
            )
        self._angles = [float(a) for a in angles]

        if isinstance(tunable_filter, str) and tunable_filter.lower() == "none":
            tunable_filter = None
        elif tunable_filter is None:
            try:
                tunable_filter = _phasing._get_tunable_filter()
            except Exception as e:
                self._logger.warning(f"Could not create the tunable filter: {e!r}")
        self._filter = tunable_filter

        self._tnfringes = tnfringes
        self._filter_mode: str | None = _read_filter_mode(tunable_filter)
        self._last_measure_tn: str | None = None
        self._spls = [
            SPL(
                camera=cam if cam is not None else "none",
                tunable_filter=tunable_filter if tunable_filter is not None else "none",
                tnfringes=tnfringes,
            )
            for cam in cams
        ]

    @property
    def spls(self) -> list[SPL]:
        """The single-camera SPL instances, in camera order."""
        return list(self._spls)

    def set_tn_fringes(self, tnfringes: str | None):
        """
        Set the tracking number of the simulated fringes, for all cameras.

        Parameters
        ----------
        tnfringes : str | None
            The tracking number of the simulated fringes measurements.
        """
        self._tnfringes = tnfringes
        for spl in self._spls:
            spl.set_tn_fringes(tnfringes)

    def set_filter_mode(self, mode: str):
        """
        Set the shared filter bandwidth mode (once, for all cameras).

        Parameters
        ----------
        mode : str
            One of "narrow", "medium", "wide", "black". The last non-black
            mode is remembered and restored after dark acquisitions.

        Raises
        ------
        ValueError
            If the mode is not valid.
        RuntimeError
            If there is no tunable filter.
        """
        from opticalib.procedures.phasing import _FILTER_BANDWIDTH_MODE as _modes

        if mode not in _modes:
            raise ValueError(
                f"Invalid filter mode: {mode}. Must be one of {list(_modes)}"
            )
        self._require_filter().set_bandwidth_mode(_modes[mode])
        if mode != "black":
            self._filter_mode = mode
        self._logger.info(f"Set filter bandwidth mode to: {mode}")

    def acquire_dark_frame(
        self,
        exptime: float,
        nframes: int = 1
    ) -> list[_ot.ImageData]:
        """
        Acquire a dark frame on every camera, with the filter in black mode.

        Parameters
        ----------
        exptime : float
            Exposure time in [s].
        nframes : int
            Number of frames to average.

        Returns
        -------
        list[ImageData]
            Copies of the dark frames (scaled to 1 s), in camera order.

        Raises
        ------
        SPLChainError
            If one or more cameras fail (the others complete anyway).
        RuntimeError
            If there is no tunable filter.

        Notes
        -----
        The filter is switched to "black" for the acquisition and the previous
        (non-black) bandwidth mode is restored afterwards, even on failure.
        """
        cur = _read_filter_mode(self._filter)
        if cur:
            self._filter_mode = cur
        self.set_filter_mode("black")
        # a camera may have been reconnected: never trust the cached exposure
        self._invalidate_exptimes()

        def _dark(k: int):
            spl = self._spls[k]
            try:
                self._set_exptime(spl, exptime)
                dark = spl._camera.acquire_frames(nframes) / exptime
            except BaseException:
                spl._curr_exptime = None
                raise
            mask = dark.mask if hasattr(dark, "mask") else _np.zeros(dark.shape)
            dark = _fits_array(data=dark, mask=mask)
            dark.header["EXPTIME"] = exptime
            dark.header["NFRAMES"] = nframes
            dark.header["ISDARK"] = True
            dark.header["ISRAW"] = True
            spl._darkFrame1sec = dark
            return dark.copy()

        try:
            with _fut.ThreadPoolExecutor(max_workers=len(self._spls)) as pool:
                futs = {k: pool.submit(_dark, k) for k in range(len(self._spls))}
            res, errs = self._collect(futs)
        finally:
            if self._filter_mode is not None:
                self.set_filter_mode(self._filter_mode)
        if errs:
            raise _chain_error("Dark frame acquisition", errs)
        return [res[k] for k in range(len(self._spls))]

    def preview_detection(
        self,
        exptime: float | None = None,
        filter_mode: str | None = None,
        wavelength: float | None = None,
        **kwargs,
    ) -> None:
        """
        Acquire one frame per camera and show the detected PSF centroid and
        crop box, one panel per camera.

        Parameters
        ----------
        exptime : float, optional
            Exposure time in [s] for all cameras. If None, the current one.
        filter_mode : str, optional
            Filter bandwidth mode to set. If None, the current one.
        wavelength : float, optional
            Wavelength to move the filter to. If None, the current one.
        **kwargs
            Passed to ``SPL.detect_psf_centroids`` (``n_psfs`` is fixed to 1).

        Returns
        -------
        tuple
            (fig, axs) of the plot.

        Raises
        ------
        SPLChainError
            If setting the exposure or acquiring fails on one or more cameras.
            A failed detection on a camera is only logged.
        """
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle

        n = len(self._spls)
        if filter_mode is not None:
            self.set_filter_mode(filter_mode)
        if exptime is not None:
            self._run(lambda k: self._set_exptime(self._spls[k], exptime), "Set exptime")
        if wavelength is not None:
            self._require_filter().move_to(wavelength)
        imgs = self._run(lambda k: self._spls[k]._camera.acquire_frames(1), "Acquire")

        ncols = min(n, 3)
        nrows = -(-n // ncols)
        fig, axs = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows), squeeze=False)
        for k, spl in enumerate(self._spls):
            ax = axs.flat[k]
            img = imgs[k]
            ax.imshow(img, cmap="gray", origin="upper")
            try:
                kw = {**kwargs, "n_psfs": 1}
                centroids = spl.detect_psf_centroids(img, **kw)
                _, boxes = spl.crop_around_centroids(img, centroids)
                for (xc, yc), (x0, y0, w, h) in zip(centroids, boxes):
                    ax.plot(xc, yc, "r+", markersize=12, markeredgewidth=2)
                    ax.add_patch(
                        Rectangle((x0, y0), w, h, fill=False, edgecolor="cyan", linewidth=1.5)
                    )
            except Exception as e:
                self._logger.warning(f"Detection failed on cam{k}: {e!r}")
            ax.set_title(f"cam{k} ({self._names[k]})")
        for ax in list(axs.flat)[n:]:
            ax.axis("off")
        plt.tight_layout()
        plt.show()
        return fig, axs

    def acquire(
        self,
        exptime: float | list[float],
        filter_mode: str | None = None,
        lambda_vector: _ot.ArrayLike | None = None,
        nframes: int = 1,
        masks: list | None = None,
    ) -> str:
        """
        Acquire an SPL sweep on all cameras at the same wavelengths.

        At each wavelength, exposure updates run in parallel with the filter
        move, then all cameras expose simultaneously, and frames are written
        in background while the sweep advances.

        Parameters
        ----------
        exptime : float | list[float]
            Base exposure time in [s]: one value for all cameras, or one per
            camera.
        filter_mode : str, optional
            Sweep filter mode ("narrow", "medium", "wide"). If None, the mode
            the filter currently has is used; if it cannot be read (or is
            black), the last mode known (set with `set_filter_mode` or read
            earlier) is used. Raises if none is known.
        lambda_vector : ArrayLike, optional
            Wavelengths between 400 and 720 nm. If None, the default vector.
        nframes : int
            Frames averaged per wavelength.
        masks : list, optional
            One mask per camera.

        Returns
        -------
        tn : str
            Tracking number of the measurement, in the SPL data folder. Data
            are in ``<tn>/cam{k}/``.

        Raises
        ------
        ValueError
            On invalid arguments, or if no sweep filter mode is known.
        SPLChainError
            If one or more cameras fail, or frame writing fails (the sweep stops
            at that wavelength). Its ``tn`` attribute is the tracking number of
            the partial measurement (None if it failed before the data folder
            was created). The filter is parked at 600 nm in any case; if that
            fails too, the failure is attached to the error as a note.
        RuntimeError
            If the sweep succeeded but the filter could not be parked at
            600 nm (the data are saved, the message gives the tracking number).
        """
        # Arguments checks
        n = len(self._spls)
        lv = SPL._validate_lambda(lambda_vector)
        gains = SPL._expgain(lv)
        if _np.ndim(exptime) == 0:
            exptimes = [float(exptime)] * n
        else:
            exptimes = [float(t) for t in exptime]
            if len(exptimes) != n:
                raise ValueError(f"Got {len(exptimes)} exposure times for {n} cameras")
        if masks is not None and len(masks) != n:
            raise ValueError(f"Got {len(masks)} masks for {n} cameras")
        if filter_mode is not None:
            mode = filter_mode
        else:
            mode = _read_filter_mode(self._filter) or self._filter_mode
        if mode is None:
            raise ValueError(
                "No filter mode for the sweep: pass `filter_mode` or call `set_filter_mode`"
            )
        if mode == "black":
            raise ValueError("Cannot sweep in black filter mode")
        filt = self._require_filter()

        # a camera may have been reconnected: never trust the cached exposure
        self._invalidate_exptimes()
        if any(spl._darkFrame1sec is None for spl in self._spls):
            self.acquire_dark_frame(exptime=1.0, nframes=1)
        self.set_filter_mode(mode)

        datapath, tn = _osu.create_data_folder(_fn.SPL_DATA_ROOT_FOLDER, True)
        print(tn)
        _osu.save_fits(_os.path.join(datapath, "lambda_vector.fits"), lv)
        cam_dirs = []
        for k, spl in enumerate(self._spls):
            d = _os.path.join(datapath, f"cam{k}")
            _os.makedirs(d, exist_ok=True)
            _osu.save_fits(_os.path.join(d, "lambda_vector.fits"), lv)
            _osu.save_fits(_os.path.join(d, "darkFrame.fits"), spl._darkFrame1sec)
            cam_dirs.append(d)
        self._logger.info(f"Starting SPLChain acquisition with tracking number: {tn}")

        cam_pool = _fut.ThreadPoolExecutor(max_workers=n)
        writer_pool = _fut.ThreadPoolExecutor(max_workers=n)
        write_futs: dict[tuple[int, float], _fut.Future] = {}
        primary: BaseException | None = None
        try:
            for i, wl in enumerate(lv):
                g = float(gains[i])
                targets = [t * g for t in exptimes]

                # exposure updates overlap the filter move (main thread)
                efuts = {
                    k: cam_pool.submit(self._set_exptime, self._spls[k], targets[k])
                    for k in range(n)
                }
                move_err = None
                try:
                    filt.move_to(wl)
                except BaseException as e:
                    move_err = e
                _, errs = self._collect(efuts)
                self._invalidate_exptimes(errs)
                if move_err is not None:
                    raise move_err
                if errs:
                    raise _chain_error(f"Exposure update at {wl} nm", errs)

                # barrier passed: all cameras expose at the same wavelength
                self._logger.info(f"Acquiring images at {wl:.1f} [nm]")
                afuts = {
                    k: cam_pool.submit(self._spls[k]._camera.acquire_frames, nframes)
                    for k in range(n)
                }
                imgs, errs = self._collect(afuts)
                self._invalidate_exptimes(errs)
                for k, img in imgs.items():
                    spl = self._spls[k]
                    write_futs[(k, wl)] = writer_pool.submit(
                        spl._save_frame,
                        cam_dirs[k],
                        img,
                        wl,
                        exptime=spl._curr_exptime,
                        expgain=g,
                        mask=masks[k] if masks is not None else None,
                        header={"CAMNAME": self._names[k]},
                    )
                if errs:
                    raise _chain_error(f"Acquisition at {wl} nm", errs)

                # do not wait for the end of the sweep to find out that the
                # disk is full: look at the writes that already finished
                werrs = self._write_errors(write_futs, only_done=True)
                if werrs:
                    raise _chain_error(f"Frame writing (up to {wl} nm)", werrs)
        except BaseException as e:
            primary = e
        finally:
            cam_pool.shutdown(wait=True)
            writer_pool.shutdown(wait=True)

        werrs = self._write_errors(write_futs)
        if primary is None and not werrs:
            self._logger.info(f"Saved tracking number: {tn}")
            self._last_measure_tn = tn
            for spl in self._spls:
                spl._last_measure_tn = tn

        park_err: Exception | None = None
        try:
            filt.move_to(600)
        except Exception as e:
            park_err = e
            self._logger.error(f"Could not park the filter at 600 nm: {e!r}")

        if primary is None and not werrs:
            if park_err is not None:
                raise RuntimeError(
                    f"The filter could not be parked at 600 nm ({park_err!r}); "
                    f"the data are saved under tracking number {tn}"
                ) from park_err
            return tn

        # failure: report everything, never drop a camera/writer/park error
        if primary is None:
            primary = _chain_error("Frame writing", werrs)
        elif werrs:
            self._logger.error(f"Frame writing failed: {werrs}")
            if isinstance(primary, SPLChainError):
                for k, e in werrs.items():
                    if primary.errors.get(k) is e:
                        continue
                    if k in primary.errors:
                        _add_note(primary.errors[k], f"frame writing also failed: {e!r}")
                    else:
                        primary.errors[k] = e
            else:
                _add_note(primary, f"frame writing failed: {werrs}")
        if park_err is not None:
            _add_note(primary, f"the filter could not be parked at 600 nm: {park_err!r}")
        if isinstance(primary, SPLChainError):
            primary.tn = tn
        raise primary

    def analysis(
        self,
        tn: str | None = None,
        max_workers: int | None = None,
        **process_kwargs,
    ) -> list[float]:
        """
        Analyse a chain measurement, one camera (one PSF) per worker.

        Parameters
        ----------
        tn : str, optional
            Tracking number. If None, the last acquired one.
        max_workers : int, optional
            Worker processes. None: one per camera; 1: sequential, in-process
            (no ``__main__`` guard needed).
        **process_kwargs
            Passed to the SPL analysis.

        Returns
        -------
        list[float]
            One piston per camera, in camera order.

        Raises
        ------
        ValueError
            If no tracking number is available, or `tnfringes` is not set.
        SPLChainError
            If the analysis of one or more cameras fails (all the cameras are
            analysed anyway); ``errors`` is keyed by camera index.

        Notes
        -----
        With ``max_workers > 1`` the workers are *spawned* processes, which
        re-import the main module: when called from a script, the call must be
        protected by the usual ``if __name__ == "__main__":`` guard.
        """
        tn = self._resolve_tn(tn)
        if self._tnfringes is None:
            raise ValueError("`tnfringes` is not set: use `set_tn_fringes` first")
        n = len(self._spls)
        paths = [self._cam_path(tn, k) for k in range(n)]
        kw = lambda k: dict(process_kwargs, n_psfs=1, angles=[self._angles[k]])
        workers = n if max_workers is None else max_workers
        res, errs = {}, {}
        if workers == 1:
            for k in range(n):
                try:
                    res[k] = _phasing._analyze_datapath(
                        paths[k], self._tnfringes, **kw(k)
                    )
                except Exception as e:
                    errs[k] = e
        else:
            with _fut.ProcessPoolExecutor(
                workers, mp_context=_mp.get_context("spawn")
            ) as pool:
                futs = {
                    k: pool.submit(
                        _phasing._analyze_datapath, paths[k], self._tnfringes, **kw(k)
                    )
                    for k in range(n)
                }
                res, errs = self._collect(futs)
        if errs:
            raise _chain_error("Analysis", errs)
        return [float(res[k][0]) for k in range(n)]

    def plot_comparison(self, tn: str | None = None, cam: int | str = "all") -> None:
        """
        Plot the template comparison for one or all cameras.

        Parameters
        ----------
        tn : str, optional
            Tracking number. If None, the last acquired one.
        cam : int | str
            Camera index, or "all".

        Raises
        ------
        ValueError
            If no tracking number is available.
        """
        tn = self._resolve_tn(tn)
        ks = range(len(self._spls)) if cam == "all" else [int(cam)]
        for k in ks:
            self._spls[k].plot_comparison(datapath=self._cam_path(tn, k), psf_n=0)

    def close(self) -> None:
        """Close the cameras (errors are logged and swallowed)."""
        for k, spl in enumerate(self._spls):
            cam = spl._camera
            if cam is not None and hasattr(cam, "close"):
                try:
                    cam.close()
                except Exception as e:
                    self._logger.error(f"Error closing cam{k}: {e!r}")

    def __enter__(self) -> "SPLChain":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    # ------------------------------------------------------------------ helpers

    def _require_filter(self):
        if self._filter is None:
            raise RuntimeError("No tunable filter available")
        return self._filter

    def _resolve_tn(self, tn: str | None) -> str:
        if tn is None:
            tn = self._last_measure_tn
        if tn is None:
            raise ValueError("No tracking number given and no measurement acquired yet")
        return tn

    @staticmethod
    def _cam_path(tn: str, k: int) -> str:
        return _os.path.join(_fn.SPL_DATA_ROOT_FOLDER, tn, f"cam{k}")

    def _invalidate_exptimes(self, which=None) -> None:
        """Forget the cached exposure (all cameras, or those in `which`)."""
        for k, spl in enumerate(self._spls):
            if which is None or k in which:
                spl._curr_exptime = None

    @staticmethod
    def _write_errors(write_futs: dict, only_done: bool = False) -> dict:
        """First write error per camera (only among finished futures if asked)."""
        errs: dict[int, BaseException] = {}
        for (k, _wl), f in list(write_futs.items()):
            if only_done and not f.done():
                continue
            e = f.exception()
            if e is not None:
                errs.setdefault(k, e)
        return errs

    @staticmethod
    def _set_exptime(spl: SPL, exptime: float) -> None:
        if spl._curr_exptime != exptime:
            spl.set_exptime(exptime)

    @staticmethod
    def _collect(futs: dict[int, _fut.Future]) -> tuple[dict, dict]:
        """Wait for all futures; return (results, errors) keyed by camera."""
        res, errs = {}, {}
        for k, f in futs.items():
            e = f.exception()
            if e is None:
                res[k] = f.result()
            else:
                errs[k] = e
        return res, errs

    def _run(self, fn, what: str) -> list:
        """Run `fn(k)` for every camera in parallel; raise SPLChainError on failure."""
        n = len(self._spls)
        with _fut.ThreadPoolExecutor(max_workers=n) as pool:
            futs = {k: pool.submit(fn, k) for k in range(n)}
        res, errs = self._collect(futs)
        if errs:
            raise _chain_error(what, errs)
        return [res[k] for k in range(n)]
