'''
Authors
  - C. Selmi: written in 2019
              rewritten in 2022
'''

import numpy as np
import logging
from skimage import measure
from matplotlib import pyplot as plt


class ROI():
    """
    Class to be used for extracting regions of interest from the image.

    HOW TO USE IT::

        from m4.utils.roi import ROI
        roi = ROI()
    """

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('ROI:')


    def roiGenerator(self, ima):
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
        self._logger.debug('Creation of roi list')
        labels = measure.label(np.invert(ima.mask))
        #from scipy import ndimage as ndi
        #labels = ndi.label(np.invert(ima.mask))[0]
        #import skimage.morphology as skm
        #pro= skm.watershed(ima, markers)
        roiList = []
        for i in range(1, 13):
            maski = np.zeros(labels.shape, dtype=bool)
            maski[np.where(labels == i)] = 1
            final_roi = np.ma.mask_or(np.invert(maski), ima.mask)
            roiList.append(final_roi)
        return roiList

    def _plotTest(self, roiList):
        plt.figure(figsize=(16, 10))
        for i in range(0,4):
            plt.subplot(2, 2, i+1)
            plt.imshow(roiList[i])
            plt.title('roiList[%d]' %i)


    def automatical_roi_selection(self, image, segment_view, ref_mirror_in):
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
        roiList = self.roiGenerator(image)

        if segment_view is True:
            if ref_mirror_in is True:
                roi_dx = roiList[1]
                roi_sx = roiList[0]
                roi_c = roiList[2]
                roi_rm = roiList[3]
                return roi_dx, roi_sx, roi_c, roi_rm
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
                return segRoiList, roi_rm
            elif ref_mirror_in is False:
                roi_seg0 = roiList[0]
                roi_seg1 = roiList[1]
                roi_seg2 = roiList[3]
                roi_seg3 = roiList[5]
                roi_seg4 = roiList[4]
                roi_seg5 = roiList[2]
                segRoiList = [roi_seg0, roi_seg1, roi_seg2, roi_seg3, roi_seg4, roi_seg5]
                return segRoiList
        