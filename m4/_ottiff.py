import os as _os
import numpy as _np
from opticalib import typings as ot
from opticalib.core import read_config as _rif
from opticalib.core.root import folders as _fn
from opticalib.ground.osutils import (
    newtn as _ts, save_fits as _sf, load_fits as _lf
)
from opticalib.ground import roi as _roi
from opticalib.dmutils import (
    iff_acquisition_preparation as _ifa, iff_processing as ipf
)


def stack_cubes(tnlist: list[str], tn_active_roi: list[int]) -> ot.CubeData:
    """
    This function stacks multiple IFF cubes into a single one, for easier processing.
    It uses the opticalib.iff_processing.stackCubes function.

    Parameters
    ----------
    tnlist: list[str]
        List of tracking numbers of the cubes to be stacked.
    tn_active_roi: list[int]
        Ordered list of the TN's active ROI, labeled as integers (e.g. [0,1] means
        TN1 has ROI 0 active, TN2 has ROI 1 active, for a total of 2 ROIs).

    Returns
    -------
    cube : CubeData
        The output cube.
    """
    if not len(tnlist) == len(tn_active_roi):
        raise ValueError("`tnlist` and `tn_active_roi` must have the same length")
    for tn, roi in zip(tnlist, tn_active_roi):
        cube = _lf(_os.path.join(_fn.INTMAT_ROOT_FOLDER, tn, "IMCube.fits"))
        _roi.roiGenerator(cube[:,:,0], len(tn_active_roi))
        ncube = []
        for img in cube.transpose(2,0,1):
            img = img*roi - _np.mean(img*roi)
            ncube.append(img)
        ncube = _np.ma.dstack(ncube)
        _sf(_os.path.join(_fn.INTMAT_ROOT_FOLDER, tn, "IMCube_roi.fits"), ncube, overwrite=True)
    final_cube = ipf.stackCubes(tnlist)
    return final_cube



class OttIffAcquisition:

    def __init__(self, dm: ot.DeformableMirrorDevice, interf: ot.InterferometerDevice):
        """The constructor"""
        if ot.isinstance_(dm, "DeformableMirrorDevice"):
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
        if ot.isinstance_(interf, "InterferometerDevice"):
            self._interf = interf
        else:
            if all(
                hasattr(interf, method)
                for method in ["acquire_map", "capture", "produce"]
            ):
                self._interf = interf
            else:
                raise TypeError(
                    """
`interf` must be an InterferometerDevice instance, or at least implements the methods:
 - acquire_map
 - capture
 - produce"""
                )


    def iffDataAcquisition(
        self,
        modesList: ot.Optional[ot.ArrayLike] = None,
        amplitude: ot.Optional[float | ot.ArrayLike] = None,
        template: ot.Optional[ot.ArrayLike] = None,
        shuffle: bool = False,
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
        self.ifc = _ifa.IFFCapturePreparation(self._dm)
        tch = self.ifc.createTimedCmdHistory(modesList, amplitude, template, shuffle)
        tn = _ts()
        if all([len(tch) == 222, self._dm._name == 'AdOpticaDm']): # TODO: make this more robust
            tchs = [tch[:111], tch[111:]] # TODO: make this more robust
            return self._newIffDataAcquisition(
                tn, tchs, modesList, amplitude, template, shuffle
            )
        else:
            self._fileAcqManager(tn, modesList, amplitude, template)
            return self._oldIffDataAcquisition(tn, tch)


    def _prepareSegmentsRoi(self):
        """
        This function prepares and labels the ROIs for the M4 (DP) segments
        """
        from opticalib.ground import roi
        
        img = self._interf.acquire_map()
        rois = roi.roiGenerator(img, 3)
        ... # TODO: complete


    def process(self, **kwargs: dict[str,ot.Any]):
        """
        Standard IFF Data processing, using the opticalib.iff_processing module.
        """
        tn = kwargs.get('tn')
        tn = [tn] if isinstance(tn, str) else tn
        for tn in enumerate(tn):
            ipf.iffProcessing(tn = tn, **kwargs)


    def _newIffDataAcquisition(
        self,
        tn: str,
        tch: ot.MatrixLike,
        modesList: ot.ArrayLike,
        amplitude: float | ot.ArrayLike,
        template: ot.ArrayLike,
    ) -> list[str]:
        """
        Modified version of the IFF data acquisition, for the particuylar case
        of the DP, which has 2 separate shells to command. The acquisition is 
        split in 2 partsa, returning a TN for each shell.
        """
        tns = []
        for k in range(2):
            stn = tn + f"_{k}"
            tns.append(stn)
            self._fileAcqManager(stn, modesList, amplitude, template)
            self._dm.uploadCmdHistory(tch[k])
            self._dm.runCmdHistory(self._interf, save=stn)
        return tns


    def _oldIffDataAcquisition(self, tn: str, tch: ot.MatrixLike) -> str:
        """
        The original IFF data acquisition function, slightly modified to fit
        in the OttIffAcquisition class.
        """
        self._dm.uploadCmdHistory(tch)
        self._dm.runCmdHistory(self._interf, save=tn)
        return tn


    def _fileAcqManager(
        self,
        tn: str,
        modesList: ot.ArrayLike,
        amplitude: float | ot.ArrayLike,
        template: ot.ArrayLike,
    ):
        """
        This function manages the creation of the folder for the IFF data,
        and the saving of the configuration parameters used for the acquisition.
        """
        info = self.ifc.getInfoToSave()
        iffpath = _os.path.join(_fn.IFFUNCTIONS_ROOT_FOLDER, tn)
        if not _os.path.exists(iffpath):
            _os.mkdir(iffpath)
        try:
            for key, value in info.items():
                if not isinstance(value, _np.ndarray):
                    tvalue = _np.array(value)
                else:
                    tvalue = value
                if key == "shuffle":
                    with open(_os.path.join(iffpath, f"{key}.dat"), "w") as f:
                        f.write(str(value))
                else:
                    _sf(_os.path.join(iffpath, f"{key}.fits"), tvalue, overwrite=True)
        except KeyError as e:
            print(f"KeyError: {key}, {e}")
        _rif.copyIffConfigFile(tn)
        for param, value in zip(
            ["modeid", "modeamp", "template"], [modesList, amplitude, template]
        ):
            if value is not None:
                _rif.updateIffConfig(tn, param, value)
