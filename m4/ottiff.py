import os as _os
import time as _time
import shutil as _sh
import numpy as _np
from opticalib import typings as ot
from opticalib.core import read_config as _rif
from opticalib.core.root import folders as _fn
from opticalib.ground import (
    osutils as _osu,
    roi as _roi,
)
from opticalib.ground.modal_decomposer import ZernikeFitter as _ZF
from opticalib.dmutils import iff_processing as ifp, iff_module as ifm
from opticalib.analyzer import cubeRebinner
from opticalib.alignment import _sc

_lf = _osu.load_fits
_sf = _osu.save_fits


class OttIffAcquisition:

    def __init__(
        self,
        dm: ot.DeformableMirrorDevice = None,
        interf: ot.InterferometerDevice = None,
    ):
        """The constructor"""
        if dm is None:
            self._dm = None
        elif ot.isinstance_(dm, "DeformableMirrorDevice"):
            self._dm = dm
        else:
            if all(
                hasattr(dm, method)
                for method in [
                    "set_shape",
                    "get_shape",
                    "uploadCmdHistory",
                    "runCmdHistory",
                ]
            ):
                self._dm = dm
            else:
                raise TypeError(
                    """
`dm` must be a DeformableMirrorDevice instance, or at least implements the methods:
 - set_shape
 - get_shape
 - uploadCmdHistory
 - runCmdHistory"""
                )
        self._interf = interf
        self._cavity = (
            _osu.load_fits(_sc.fitting_surface)
            if not _sc.fitting_surface == ""
            else None
        )


    def acquireDpIff(
        self, separate_shells: bool = True, **kwargs: dict[str, ot.Any]
    ) -> str | list[str]:
        """
        This is the user-lever function for the acquisition of the IFF data, given a
        deformable mirror and an interferometer.

        This is an expansion of the `opticalib.iff_module.iffDataAcquisition` function,
        customized for the M4 DP segmented system, in which ROIs play a fundamental role.

        Except for the devices, all the arguments are optional, as, by default, the
        values are taken from the `configuration.yaml` file, from the IFFUNCTION section.

        Parameters
        ----------
        separate_shells: bool , optional
            if True, acquire the shells separately
        modesList: ArrayLike , optional
            list of modes index to be measured, relative to the command matrix to be used
        amplitude: float | ArrayLike, optional
            command amplitude
        template: ArrayLike , oprional
            template file for the command matrix
        shuffle: bool , optional
            if True, shuffle the modes before acquisition

        Returns
        -------
        tns: list[str]
            The tracking number of the dataset acquired, saved in the OPDImages folder
        """
        ifconfig = _rif.getIffConfig("IFFUNC")
        amplitude = kwargs.pop("amplitude", _np.full(222, ifconfig.get("amplitude")))
        modeslist = [_np.arange(0, 111, 1), _np.arange(111, 222, 1)]
        ampvec = [amplitude[:111], amplitude[111:]]
        if not separate_shells:
            modeslist = _np.arange(0, 222, 1)
            ampvec = amplitude
            return ifm.iffDataAcquisition(
                dm=self._dm,
                interf=self._interf,
                modesList=modeslist,
                amplitude=ampvec,
                **kwargs,
            )
        else:
            tns = []
            for k in range(2):
                tns.append(
                    ifm.iffDataAcquisition(
                        dm=self._dm,
                        interf=self._interf,
                        modesList=modeslist[k],
                        amplitude=ampvec[k],
                        **kwargs,
                    )
                )
            return tns

    # def process(self, **kwargs: dict[str, ot.Any]):
    #     """
    #     Standard IFF Data processing, using the opticalib.iff_processing module.
    #     """
    #     tn = kwargs.pop("tn", None)
    #     tn = [tn] if isinstance(tn, str) else tn
    #     for _, t in enumerate(tn):
    #         ifp.process(tn=t, **kwargs)

    def dpIffRoiProcessing(
        self,
        tns: list[str],
        activeRoiID: list[int],
        detrend: bool = False,
        roinull: bool = False,
        rebin: int = 1,
        fitting_mask = None
    ):
        """
        This function groups together some image manipulations in presence of
        Region Of Interest (ROI). We assume that we have in the image an
        activeRoi, with given I, corresponding to the region of the actuated
        segment; and one or more auxiliaryRois, corresponding to the regions
        of the non-actuated segments, to be used as an optical reference.

        Parameters
        ----------
        tn: str
            The tracking number of the dataset to be processed.
        activeRoiID: int
            The ID of the active ROI, corresponding to the actuated segment.
        detrend: bool, optional
            If True, perform a tilt detrend over the activeRoi, using the auxRoi as reference.
        roinull: bool, optional
            If True, set to zero the pixels outside the activeRoi.
        rebin: int, optional
            Rebin factor for the final cube.

        Returns
        -------
        newtn: str
            The tracking number of the new processed dataset.
        """
        from opticalib.ground.geometry import draw_circular_pupil

        #mask = draw_circular_pupil(shape=(2000,2000), radius=990)

        if isinstance(tns, str):
            tns = [tns]
        ifp.process(tn=tns, save=True, rebin=rebin)
        newtns = []
        for j, tn in enumerate(tns):
            newtn = ifp.cubeRoiProcessing(tn, activeRoiID[j], detrend, roinull, rebin, fitting_mask) #, fitting_mask=self._cavity)
            newtns.append(newtn)
            # _time.sleep(1) # WHY? 
        return newtns

    # def _tiltDetrend(
    #     self, img: ot.ImageData, activeRoi: ot.MaskData, auxRois: ot.MaskData
    # ):
    #     """
    #     This method detrends an image by fitting and subtracting PTT
    #     (Zernike modes 1, 2, 3) from auxiliary regions of interest (ROIs), and
    #     applies the correction to both the ROIs and the active ROI area.

    #     Parameters
    #     ----------
    #     img : ot.ImageData
    #         The input image data to be detrended.
    #     activeRoi : ot.MaskData
    #         Mask defining the active region of interest where mean correction
    #         will be applied.
    #     auxRois : ot.MaskData
    #         Collection of auxiliary ROI masks used for tilt estimation. Each ROI
    #         is processed independently.

    #     Returns
    #     -------
    #     detrendedImg : ot.ImageData
    #         The detrended image with tilt removed. For each auxiliary ROI, the
    #         fitted surface is subtracted from the entire image, and the mean of
    #         the fitted surface is subtracted from areas outside the active ROI.

    #     Notes
    #     -----
    #     The method uses Zernike polynomials [1, 2, 3] which typically correspond
    #     to piston, tip, and tilt modes for the tilt removal operation.
    #     """
    #     detrendedImg = img.copy()
    #     for roi in auxRois:
    #         r = roi.copy()
    #         r2rImage = _np.ma.masked_array(img.data, mask=r)
    #         surf2Remove = self._zern.makeSurface([1, 2, 3], r2rImage)
    #         s2rmean = surf2Remove.mean()
    #         surf2Remove.data[activeRoi == 0] = 0
    #         surf2Remove.mask[activeRoi == 0] = False

    #         detrendedImg -= surf2Remove
    #         detrendedImg -= s2rmean
    #     return detrendedImg

    # def roizern2(self, img, z2fit, auxmask =None, roiid=None, local =True, roiimg = None):
    #     if ((roiid is not None) and (roiimg is None)):  #
    #         roiimg = _roi.roiGenerator(img) #non è disponibile un parametro passato per dire QUANTE roi cercare. funziona anche senza?
    #         nroi = len(roiid)
    #     else:
    #         nroi=1
    #     if auxmask is None:
    #         auxmask2use = img.mask
    #     else:
    #         auxmask2use = auxmask
    #     zcoeff = _np.zeros([nroi, len(z2fit)])
    #     zsurf  = []
    #     for i in range(nroi):
    #         img2fit = _np.ma.masked_array(img.data, roiimg[i])
    #         cc, _ = _zern.zernikeFitAuxmask(img2fit, auxmask2use, z2fit)
    #         zcoeff[i,:] = cc
    #     if local is False:
    #         zcoeff = zcoeff.mean(axis=0)
    #     return zcoeff

    # def _tiltDetrend2(self, img, auxmask, roi2Calc, roi2Remove, roiimg=None):
    #     '''
    #     computes the Zernikes (PTT only)  over the roi2Calc, then produces the corresponding shape over the roi2Remove mask, and subtract it.
    #     USAGE:
    #     v=self._tiltDetrend(imgf,mm, [0],[1]) # 0 is the Id of the non-active ROI; 1 is the Id of the active ROI
    #     '''
    #     if roiimg is None:
    #         roiimg = _roi.roiGenerator(img)
    #     if auxmask is None:
    #         auxmask = img.mask
    #     else:
    #         auxmask = auxmask
    #     zcoeff = self.roizern(img, [1,2,3], auxmask, roiid=roi2Calc, local=True, roiimg=roiimg) #returns the global PTT evaluated over the roi2Calc areas
    #     am = _np.ma.masked_array(auxmask, mask=auxmask==0)
    #     _, zmat = _zern.zernikeFit(am,[1,2,3]) #returns the ZernMat created over the entire circular pupil
    #     surf2Remove = _zern.zernikeSurface(am, zcoeff[roi2Calc,:], zmat)
    #     #surf2Remove[roiimg[roi2Remove==0]] =0
    #     #surf2Remove = np.ma.masked_array(surf2remove.data, roi2Remove.mask)
    #     detrendedImg = img - surf2Remove
    #     return detrendedImg
