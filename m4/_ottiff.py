import shutil as _sh
import os as _os
import numpy as _np
from opticalib import typings as ot
from opticalib.core import read_config as _rif
from opticalib.core.root import folders as _fn
from opticalib.ground import (
    osutils as _osu,
    roi as _roi,
    zernike as _zern
)
from opticalib.dmutils import (
    iff_acquisition_preparation as _ifa, 
    iff_processing as ipf,
    iff_module as ifm
)
from opticalib.analyzer import cubeRebinner

_lf = _osu.load_fits
_sf = _osu.save_fits

def stackRoiCubes(tnlist: list[str], tn_active_roi: list[int]) -> ot.CubeData:
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
    for r, tn in enumerate(tnlist):
        cube = _lf(_os.path.join(_fn.INTMAT_ROOT_FOLDER, tn, "IMCube.fits"))
        rois = _roi.roiGenerator(cube[:,:,0], len(tn_active_roi))
        ncube = []
        rr = tn_active_roi[r]
        for img in cube.transpose(2,0,1):
            fin_img = img.copy()
            # POSSIBILE ASTRAZIONE: la shell attiva è la `r`-esima, mentre quelle
            # non attive sono tutte le altre, e quindi ciclo su ROIS con 
            # if roi==r continue, così che si crea una lista di non active shells
            # da poi unire nell'immagine finale
            nonactiveroi = rois[1-rr]
            fin_img[nonactiveroi==0] = 0
            ncube.append(fin_img)
        ncube = _np.ma.dstack(ncube)
        # fare la maschera del cubo, così che venga salvata
        # ... TODO
        import subprocess
        original_path = _os.path.join(_fn.INTMAT_ROOT_FOLDER, tn, "IMCube.fits")
        renamed_path = _os.path.join(_fn.INTMAT_ROOT_FOLDER, tn, "og_IMCube.fits")
        subprocess.run(['mv', original_path, renamed_path])
        _sf(_os.path.join(_fn.INTMAT_ROOT_FOLDER, tn, "IMCube.fits"), ncube, overwrite=True)
    final_cube = ipf.stackCubes(tnlist)
    return final_cube



class OttIffAcquisition:

    def __init__(self, dm: ot.DeformableMirrorDevice = None, interf: ot.InterferometerDevice = None):
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
        if interf is None:
            self._interf = None
        elif ot.isinstance_(interf, "InterferometerDevice"):
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
        


    def acquireDpIff(
        self,
        dp: bool = True,
        **kwargs: dict[str,ot.Any]
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
        amplitude = kwargs.pop('amplitude', _np.full(222, 500e-9))
        if dp is True:
            modeslist = [_np.arange(0,111,1), _np.arange(111,222,1)]
            ampvec = [amplitude[:111], amplitude[111:]]
            tns = []
            for k in range(2):
                tns.append(ifm.iffDataAcquisition(
                    dm=self._dm, interf=self._interf, modesList=modeslist[k], amplitude=ampvec[k], **kwargs
                ))
            return tns
        else:
            return ifm.iffDataAcquisition(dm=self._dm, interf=self._interf, amplitude=amplitude, **kwargs)

    def process(self, **kwargs: dict[str,ot.Any]):
        """
        Standard IFF Data processing, using the opticalib.iff_processing module.
        """
        tn = kwargs.pop('tn', None)
        tn = [tn] if isinstance(tn, str) else tn
        for _,t in enumerate(tn):
            ipf.process(tn = t, **kwargs)


    def iffRoiProcessing(self, 
                          tns: list[str], 
                          auxmask: ot.ImageData,
                          activeRoiID: list[int],
                          auxRoiID: list[int],
                          detrend: bool = False,
                          roinull: bool = False,
                          roicosmetic: bool = False,
                          rebin:int = 1
                          ):
        """
        This function groups together some image manipulations in presence of 
        Region Of Interest (ROI). We assume that we have in the image an
        activeRoi, with given I, corresponding to the region of the actuated
        segment; and one or more auxiliaryRosi, corresponding to the regions
        of the non-actuated segments, to be used as an optical reference.

        Parameters
        ----------
        tn: str
            The tracking number of the dataset to be processed.
        auxmask: 2D array-like
            The auxiliary mask, with value 1 in the region of interest, 0 elsewhere.
        activeRoiID: list[int]
            The ID of the active ROI, corresponding to the actuated segment.
        auxRoiID: list[int]
            The ID(s) of the auxiliary ROIs, corresponding to the non-actuated segments.
        detrend: bool, optional
            If True, perform a tilt detrend over the activeRoi, using the auxRoi as reference.
        roinull: bool, optional
            If True, set to zero the pixels outside the activeRoi.
        roicosmetic: bool, optional
            If True, perform cosmetic operations on the image (e.g. remove bad pixels).
        """
        newtn = _osu.newtn() # ??
        cube = _lf(_os.path.join(_fn.INTMAT_ROOT_FOLDER, tns, "IMCube.fits")).transpose(2,0,1)
        mat = _lf(_os.path.join(_fn.INTMAT_ROOT_FOLDER, tns, "cmdMatrix.fits"))
        modesvec = _lf(_os.path.join(_fn.INTMAT_ROOT_FOLDER, tns, "modesVector.fits"))
        nframes = cube.shape[0]
        newcube = []
        for i in range(nframes):
            rois = _roi.roiGenerator(cube[i])
            img = cube[i]
            if detrend is not False:
                v = self._tiltDetrend(img, auxmask, auxRoiID, activeRoiID, rois)
            if roinull is not False:
                v = self.roinull(v, auxRoiID, rois)
            if roicosmetic is not False:
                v = self.roicosmetics(v)
            newcube.append(v)
        newcube = _np.ma.masked_array(newcube)
        newcube = cubeRebinner(newcube.transpose(1,2,0), rebin)
        save_path = _os.path.join(_fn.INTMAT_ROOT_FOLDER, newtn)
        if not _os.path.exists(save_path):
            _os.makedirs(save_path)
        _osu.save_fits(_os.path.join(save_path, "IMCube.fits"), newcube, overwrite=True)
        _osu.save_fits(_os.path.join(save_path, "cmdMatrix.fits"), mat, overwrite=True)
        _osu.save_fits(_os.path.join(save_path, "modesVector.fits"), modesvec, overwrite=True)
        oflag = _os.path.join(_fn.INTMAT_ROOT_FOLDER, tns, "flag.txt")
        nflag = _os.path.join(save_path, "flag.txt")
        _sh.copy(oflag, nflag)
        return newtn 

    def roicosmetics(self, img, params=None):
        return img
    

    def roinull(self, img, roiId, roiimg = None):
        if roiimg is None:
            roiimg = _roi.roiGenerator(img)
        for i in roiId:
            img[roiimg[i]==0] = 0
        return img

    def roizern(self, img, z2fit, auxmask =None, roiid=None, local =True, roiimg = None):
        if roiid is None:
            nroi = 1
        else:
            if roiimg is None:
                print('roizern: searching rois')
                roiimg = _roi.roiGenerator(img) #non è disponibile un parametro passato per dire QUANTE roi cercare. funziona anche senza?
            nroi = len(roiid)

        print('nroi')
        print(nroi)
        if auxmask is None:
            auxmask2use = img.mask
        else:
            auxmask2use = auxmask
        #zcoeff = []
        zcoeff =  _np.zeros([nroi, len(z2fit)])
        zsurf = []
        for i in range(nroi):
            img2fit = _np.ma.masked_array(img.data, roiimg[i])
            cc, zmat = _zern.zernikeFitAuxmask(img2fit, auxmask2use, z2fit)
            #qui implementare il return della zsurface
            #zcoeff.append(cc)
            zcoeff[i] = cc  #was zcoeff[i,:] = cc
            zsurf.append(zmat)
        #zcoeff = _np.array(zcoeff)
        if (local is False) and (nroi >1):
            print('Option -global- selected')
            zcoeff = zcoeff.mean(axis=0)
        print('zcoeff roizern')
        print(zcoeff)
        return zcoeff

    def _tiltDetrend(self, img, auxmask, roi2Calc, roi2Remove, roiimg=None, zsurf=None):
        '''
        computes the Zernikes (PTT only)  over the roi2Calc, then produces the corresponding shape over the roi2Remove mask, and subtract it.
        USAGE:
        v=self._tiltDetrend(imgf,mm, [0],[1]) # 0 is the Id of the non-active ROI; 1 is the Id of the active ROI
        '''
        if roiimg is None:
            roiimg = _roi.roiGenerator(img)
        if auxmask is None:
            auxmask = img.mask
        else:
            auxmask = auxmask
        zcoeff = self.roizern(img, [1,2,3], auxmask, roiid=roi2Calc, local=False, roiimg=roiimg) #returns the global PTT evaluated over the roi2Calc areas
        am = _np.ma.masked_array(auxmask, mask=auxmask==0)
        _, zmat = _zern.zernikeFit(am,[1,2,3]) #returns the ZernMat created over the entire circular pupil
        #      zmat = zsurf[roi2Remove[0]]
        surf2Remove = _zern.zernikeSurface(am, zcoeff[0], zmat)  #was zcoeff[roi2Calc,:]
        #surf2Remove[roiimg[roi2Remove==0]] =0
        #surf2Remove = np.ma.masked_array(surf2remove.data, roi2Remove.mask)
        detrendedImg = img - surf2Remove
        detrendedImg[roiimg[roi2Remove[0]]==0] -= detrendedImg[roiimg[roi2Remove[0]]==0].mean()
        return detrendedImg

    def roizern2(self, img, z2fit, auxmask =None, roiid=None, local =True, roiimg = None):
        if ((roiid is not None) and (roiimg is None)):  #
            roiimg = _roi.roiGenerator(img) #non è disponibile un parametro passato per dire QUANTE roi cercare. funziona anche senza?
            nroi = len(roiid)
        else:
            nroi=1
        if auxmask is None:
            auxmask2use = img.mask
        else:
            auxmask2use = auxmask
        zcoeff = _np.zeros([nroi, len(z2fit)])
        zsurf  = []
        for i in range(nroi):
            img2fit = _np.ma.masked_array(img.data, roiimg[i])
            cc, _ = _zern.zernikeFitAuxmask(img2fit, auxmask2use, z2fit)
            zcoeff[i,:] = cc
        if local is False:
            zcoeff = zcoeff.mean(axis=0)
        return zcoeff

    def _tiltDetrend2(self, img, auxmask, roi2Calc, roi2Remove, roiimg=None):
        '''
        computes the Zernikes (PTT only)  over the roi2Calc, then produces the corresponding shape over the roi2Remove mask, and subtract it.
        USAGE:
        v=self._tiltDetrend(imgf,mm, [0],[1]) # 0 is the Id of the non-active ROI; 1 is the Id of the active ROI
        '''
        if roiimg is None:
            roiimg = _roi.roiGenerator(img)
        if auxmask is None:
            auxmask = img.mask
        else:
            auxmask = auxmask
        zcoeff = self.roizern(img, [1,2,3], auxmask, roiid=roi2Calc, local=True, roiimg=roiimg) #returns the global PTT evaluated over the roi2Calc areas
        am = _np.ma.masked_array(auxmask, mask=auxmask==0)
        _, zmat = _zern.zernikeFit(am,[1,2,3]) #returns the ZernMat created over the entire circular pupil
        surf2Remove = _zern.zernikeSurface(am, zcoeff[roi2Calc,:], zmat)
        #surf2Remove[roiimg[roi2Remove==0]] =0
        #surf2Remove = np.ma.masked_array(surf2remove.data, roi2Remove.mask)
        detrendedImg = img - surf2Remove
        return detrendedImg


