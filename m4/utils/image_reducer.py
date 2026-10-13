import logging
import numpy as np
from skimage import measure
from opticalib.ground import modal_decomposer as _md

# ---------------------------------------------------------------------------
# Legacy ROI helpers restored verbatim from the removed m4.utils.roi
# (roiGenerator, automatical_roi_selection).
# ---------------------------------------------------------------------------


def _roiGenerator(ima):
    '''
    Parameters
    ----------
        ima: numpy masked array
            image

    Returns
    -------
        roiList: list
            list of the first 12 roi found in the image

    .. note::

        roiList[3] = RM roi for alignement, roiList[3] = central roi for segment

    '''
    labels = measure.label(np.invert(ima.mask))
    roiList = []
    for i in range(1, 13):
        maski = np.zeros(labels.shape, dtype=bool)
        maski[np.where(labels == i)] = 1
        final_roi = np.ma.mask_or(np.invert(maski), ima.mask)
        roiList.append(final_roi)
    return roiList

def _automatical_roi_selection(image, segment_view, ref_mirror_in):
    '''
    Parameters
    ----------
    image: numpy masked array
        image to be analyzed
    segment_view = boolean
        in the ott is in segment view configuration it is True,
        else False
    RM_in = boolean
        if reference mirror is inside the image it is True,
        else False
    '''
    roiList = _roiGenerator(image)

    if segment_view is True:
        if ref_mirror_in is True:
            roi_dx = roiList[1]
            roi_sx = roiList[0]
            roi_c = roiList[2]
            roi_rm = roiList[3]
        elif ref_mirror_in is False:
            roi_dx = roiList[2]
            roi_sx = roiList[1]
            roi_c = roiList[3]
            roi_rm = roiList[0]
        return roi_dx, roi_sx, roi_c, roi_rm

    elif segment_view is False:
        if ref_mirror_in is True:
            roi_seg0 = roiList[0]
            roi_seg1 = roiList[1]
            roi_seg2 = roiList[3]
            roi_seg3 = roiList[6]
            roi_seg4 = roiList[5]
            roi_seg5 = roiList[2]
            segRoiList = [roi_seg0, roi_seg1, roi_seg2, roi_seg3, roi_seg4, roi_seg5]
            roi_rm = roiList[4]
        elif ref_mirror_in is False:
            roi_seg0 = roiList[0]
            roi_seg1 = roiList[1]
            roi_seg2 = roiList[3]
            roi_seg3 = roiList[5]
            roi_seg4 = roiList[4]
            roi_seg5 = roiList[2]
            segRoiList = [roi_seg0, roi_seg1, roi_seg2, roi_seg3, roi_seg4, roi_seg5]
            roi_rm = None
        return segRoiList, roi_rm


class TipTiltDetrend:
    """
    Class for removal of image's tip and tilt

    HOW TO USE IT::

        from m4.utils.image_reducer import TipTiltDetrend
        TT = TipTiltDetrend()
    """

    def __init__(self):
        """The constructor"""
        self._logger = logging.getLogger("TIP_TILT_DETREND:")

    def segment_view_tiptilt_detrend(self, image):
        """
        Parameters
        ----------
            image: numpy masked array
                    image to be analyzed

        Returns
        -------
                image_ttr: numpy array
                         image without tip and tilt on central segment
        """
        self._logger.debug("Removal of tip-tilt from central segment")
        roi_dx, roi_sx, roi_c, roi_rm = _automatical_roi_selection(
            image, True, False
        )

        dx_image = np.ma.masked_array(image.data, mask=roi_dx)
        sx_image = np.ma.masked_array(image.data, mask=roi_sx)
        coef_dx, _ = _md.ZernikeFitter(dx_image).fit(dx_image, np.arange(10) + 1)
        coef_sx, _ = _md.ZernikeFitter(sx_image).fit(sx_image, np.arange(10) + 1)
        coef = np.array(
            [(coef_dx[0] + coef_sx[0]) / 2.0, (coef_dx[1] + coef_sx[1]) / 2.0]
        )

        central_image = np.ma.masked_array(image.data, mask=roi_c)
        zf_c = _md.ZernikeFitter(central_image)
        _, mat = zf_c.fit(central_image, [2, 3])
        ima_ttr = zf_c.make_surface([2, 3], central_image, coeffs=coef, mat=mat)
        return ima_ttr

    def central_view_tiptilt_detrend(self, image, segment_ind):
        """
        Parameters
        ----------
            image: numpy masked array
                    image to be analyzed
            segment_ind: int
                    segment number to be analyzed

        Returns
        -------
                image_ttr: numpy array
                         image without tip and tilt on segment ind
        """
        segRoiList, roi_rm = _automatical_roi_selection(image, False, True)
        final_roi = segRoiList[segment_ind]
        final_image = np.ma.masked_array(image, mask=final_roi)
        del segRoiList[segment_ind]

        tip_list = []
        tilt_list = []
        for roi in segRoiList:
            ima = np.ma.masked_array(image, mask=roi)
            cc, _ = _md.ZernikeFitter(ima).fit(ima, np.arange(10) + 1)
            tip_list.append(cc[0])
            tilt_list.append(cc[1])
        tip = np.mean(np.array(tip_list))
        tilt = np.mean(np.array(tilt_list))
        coef = np.array([tip, tilt])

        zf_f = _md.ZernikeFitter(final_image)
        _, mat = zf_f.fit(final_image, [2, 3])
        ima_ttr = zf_f.make_surface([2, 3], final_image, coeffs=coef, mat=mat)
        return ima_ttr
