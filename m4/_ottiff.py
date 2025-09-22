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
    iff_acquisition_preparation as _ifa, 
    iff_processing as ipf,
    iff_module as ifm
)


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
#            active_shell = _np.ma.masked_array(uimg.data, mask=rois[r])
#            nonactive_shell = _np.ma.masked_array(uimg.data, mask=rois[1-r])
#            active_shell -= active_shell.mean()
#            nonactive_shell[rois[1-r]==0] = 0
#            total_data = _np.zeros(uimg.shape)
#            total_data[rois[1-r]] = nonactive_shell[rois[1-r]]
#            total_data[rois[r]] = active_shell[rois[r]]
#            total_mask = _np.logical_xor(rois[r], rois[1-r])
#            fin_img = _np.ma.masked_array(total_data, mask=_np.invert(total_mask))
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
        if dp is True:
            modeslist = [_np.arange(0,111,1), _np.arange(111,222,1)]
            ampvec = [amplitude[:111], amplitude[111:]]
            tns = []
            for k in range(2):
                tns.append(ifm.iffDataAcquisition(
                    dm=self._dm, interf=self._interf, modesList=modeslist[k], **kwargs
                ))
            return tns
        else:
            return ifm.iffDataAcquisition(dm=self._dm, interf=self._interf, **kwargs)

    def process(self, **kwargs: dict[str,ot.Any]):
        """
        Standard IFF Data processing, using the opticalib.iff_processing module.
        """
        tn = kwargs.pop('tn', None)
        tn = [tn] if isinstance(tn, str) else tn
        for _,t in enumerate(tn):
            ipf.process(tn = t, **kwargs)


    def _prepareSegmentsRoi(self):
        """
        This function prepares and labels the ROIs for the M4 (DP) segments
        """
        from opticalib.ground import roi
        
        img = self._interf.acquire_map()
        rois = roi.roiGenerator(img, 3)
        ... # TODO: complete


    def iff_roiProcessing(tn, auxmask, activeRoiID, auxRoiID, detrend=False, roinull=False, roicosmetic=False):
        '''
        This function groups together some image manipulations in presence of Region Of Interest (ROI). We assume thatwe have in the image an activeRoi, with given I, corresponding to the region of the actuated segment; and one or more auxiliaryRosi, corresponding to the regions of the non-actuated segments, to be used as an optical reference.

        Parameters
        ----------
        tn     :    str
            The tracking number of the dataset (expected: the IFF cube, saved by the iff_process script)
        auxmask       :  
            The circular mask to create the Zernike modes
        activeRoiID   :
            The ID of the region 
        auxRoiID      :
            
        detrend       :
            
        roinull       :
            
        roicosmetic   :
            
        Returns
        -------

        '''
        #qui legge il cubo da TN, trova il num di immagini e altro.
        # supponiamo che il cubo si chiami cube e che possa indicizzarlo come cube[i] e che contenga nframes.
        # supponiamo tutto il cubo sia mascherato con la maschera intersezione. la master-mask corrisponde quindi alla maschera di un frame qualsiasi
        #      master_mask = TBD
        nrois2search = len([activeRoiID,auxRoiID])
        rois = roi.roiGenerator(master_mask,nrois2search)
        newcube = []
        for i in range(nframes):
            img = cube[i]
            if detrend is not False:
                v=tiltDetrend(img,auxmask, auxRoiID,activeRoiID, rois)
            if roinull is not False:
                v = roinull(v, auxRoiID, rois)
            if roicosmetic is not False:
                v = roicosmetics(img)
            newcube.append(v)
            

        newcube = np.ma.masked_array(newcube)
        #qui salvare il newcube
                    

    def roicosmetics(img, params=None):

        return img

    def roiId(img):
        roiimg = _roi.roiGenerator(img)
        print('N. of ROIs found:')
        print(len(roiimg)
            for i in roiimg:
                imshow(i)
                title(str(i))


    def roinull(img, roi2null, roiimg = None):
        if roiimg is None:
            roiimg = _roi.roiGenerator(img)
        for i in roi2null:
            img[i==0]=0
        return img

    def roizern(img, z2fit, auxmask =None, roiid=None, local =True, roiimg = None):
        if ((roiid is not None) and (roiimg is None)):  #
            roiimg = _roi.roiGenerator(img) #non è disponibile un parametro passato per dire QUANTE roi cercare. funziona anche senza?
            nroi = len(roiid)
        else:
            nroi=1
        if auxmask is None:
            auxmask2use = img.mask
        else:
            auxmask2use = auxmask
        zcoeff = np.zeros([nroi, len(z2fit)])
        zsurf  = []
        for i in range(nroi):
            img2fit = np.ma.masked_array(img.data, roiimg[i])
            cc, _ =zern.zernikeFitAuxmask(img2fit, auxmask2use, z2fit)
            zcoeff[i,:] = cc
        if local is False:
            zcoeff = zcoeff.mean(axis=0)
        return zcoeff

    def _tiltDetrend(img, auxmask, roi2Calc, roi2Remove, roiimg=None):
       '''
        computes the Zernikes (PTT only)  over the roi2Calc, then produces the corresponding shape over the roi2Remove mask, and subtract it.
        USAGE:
        v=tiltDetrend(imgf,mm, [0],[1]) # 0 is the Id of the non-active ROI; 1 is the Id of the active ROI
        '''
        if roiimg is None:
            roiimg = _roi.roiGenerator(img)
        zcoeff = roizern(img, [1,2,3], auxmask, roiid=roi2Calc, local=True, roiimg) #returns the global PTT evaluated over the roi2Calc areas
        am = np.ma.masked_array(auxmask, auxmask==0)
        _, zmat = zern.zernikeFit(am,[1,2,3]) #returns the ZernMat created over the entire circular pupil
        surf2Remove = zern.zernikeSurface( am,zcoeff[roi2Calc,:], zmat)
        #surf2Remove[roiimg[roi2Remove==0]] =0
        #surf2Remove = np.ma.masked_array(surf2remove.data, roi2Remove.mask)
        detrendedImg = imgf - surf2Remove
        return detrendedImg


