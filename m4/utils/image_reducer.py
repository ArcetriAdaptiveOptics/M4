import logging
import numpy as np
from opticalib.ground import roi, zernike


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
        self.roi = roi

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
        roi_dx, roi_sx, roi_c, roi_rm = self.roi.automatical_roi_selection(
            image, True, False
        )

        dx_image = np.ma.masked_array(image.data, mask=roi_dx)
        sx_image = np.ma.masked_array(image.data, mask=roi_sx)
        coef_dx, mat_dx = zernike.zernikeFit(dx_image, np.arange(10) + 1)
        coef_sx, mat_sx = zernike.zernikeFit(sx_image, np.arange(10) + 1)
        coef = np.array(
            [(coef_dx[0] + coef_sx[0]) / 2.0, (coef_dx[1] + coef_sx[1]) / 2.0]
        )

        central_image = np.ma.masked_array(image.data, mask=roi_c)
        cc, mat = zernike.zernikeFit(central_image, np.arange(10) + 1)
        ima_ttr = zernike.zernikeSurface(central_image, coef, mat[:, 1:3])
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
        segRoiList, roi_rm = self.roi.automatical_roi_selection(image, False, True)
        final_roi = segRoiList[segment_ind]
        final_image = np.ma.masked_array(image, mask=final_roi)
        del segRoiList[segment_ind]

        tip_list = []
        tilt_list = []
        for roi in segRoiList:
            ima = np.ma.masked_array(image, mask=roi)
            cc, mat = zernike.zernikeFit(ima, np.arange(10) + 1)
            tip_list.append(cc[0])
            tilt_list.append(cc[1])
        tip = np.mean(np.array(tip_list))
        tilt = np.mean(np.array(tilt_list))
        coef = np.array([tip, tilt])

        cc, mat = zernike.zernikeFit(final_image, np.arange(10) + 1)
        ima_ttr = zernike.zernikeSurface(final_image, coef, mat[:, 1:3])
        return ima_ttr
