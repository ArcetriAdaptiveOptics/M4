"""
Author(s)
    - Chiara Selmi: written in 2019
                    rewritten in 2022
    - Pietro Ferraiuolo: modified in 2024
"""
import logging
import numpy as np
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
#       from scipy import ndimage as ndi
#       labels = ndi.label(np.invert(ima.mask))[0]
#       import skimage.morphology as skm
#       pro= skm.watershed(ima, markers)
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

    def single_segment_mask(self, image, apply:bool=True):
        """
        Given an interferometer image of a selected, non-masked, M4's segment,
        automatically masks the image so that only the pointed segment is visible.

        Parameters
        ----------
        image : masked ndarray
            Interferometer images of a segment, in which adjacent segments are
            visible.
        apply : bool, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        img_out : masked ndarray
            Input image masked so that only the principal segment is visible.
        """
        roilist = measure.label(np.invert(image.mask))
        segments = self._find_big_segments_roi(roilist)
        max_area = max([seg.area for seg in segments])
        act_seg = segments[next(i for i, obj in enumerate(segments) \
                                                      if obj.area == max_area)]
        maski = np.zeros(roilist.shape, dtype=bool)
        maski[np.where(roilist == act_seg.label)] = 1
        final_roi = np.ma.mask_or(np.invert(maski), image.mask)
        if apply:
            out = np.ma.masked_array(data=image, mask=final_roi)
        else:
            out = final_roi
        return out

    def adjacent_segments_mask(self,image, apply:bool=True):
        """
        Given an interferometer image of a selected, non-masked, M4's segment,
        automatically masks the image so that only the pointed segment and it's
        two adjacent ones are visible.

        Parameters
        ----------
        image : masked ndarray
            Interferometer images of a segment, in which all segments in frame
            are visible.
        apply : bool, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        img_out : masked ndarray
            Input image masked so that only the three principal segments are
            visible (the pointed one and its adjacents).
        """
        roilist = measure.label(np.invert(image.mask))
        segments = self._find_big_segments_roi(roilist)
        maski = mask1 = mask2 = mask3 = np.zeros(roilist.shape, dtype=bool)
        mask1[(np.where(roilist == segments[0].label))] = 1
        mask2[(np.where(roilist == segments[1].label))] = 1
        mask3[(np.where(roilist == segments[2].label))] = 1
        maski = np.ma.mask_or(np.ma.mask_or(mask1, mask2), mask3)
        final_roi = np.ma.mask_or(np.invert(maski), image.mask)
        if apply:
            out = np.ma.masked_array(data=image, mask=final_roi)
        else:
            out = final_roi
        return out

    def only_segments_mask(self, image, apply:bool=True):
        """


        Parameters
        ----------
        image : TYPE
            DESCRIPTION.
        apply : bool, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        img_out : TYPE
            DESCRIPTION.

        """
        roilist = measure.label(np.invert(image.mask))
        segments = self._find_all_segments_roi(roilist)
        masks = [np.zeros(roilist.shape, dtype=bool)]
        for i, seg in enumerate(segments):
            mask = np.zeros(roilist.shape, dtype=bool)
            mask[(np.where(roilist == seg.label))] = 1
            masks.append(mask)
        for i in range(1, len(masks)):
            masks[0] = np.ma.mask_or(masks[0], masks[i])
        final_roi = np.ma.mask_or(np.invert(masks[0]), image.mask)
        if apply:
            out = np.ma.masked_array(data=image, mask=final_roi)
        else:
            out = final_roi
        return out

    def _find_all_segments_roi(self, roilist):
        """


        Parameters
        ----------
        roilist : TYPE
            DESCRIPTION.

        Returns
        -------
        regions : TYPE
            DESCRIPTION.

        """
        regions = measure.regionprops(roilist)
        for i,region in enumerate(regions):
            if region.area < 10000: #TODO - Find good criteria
                regions.pop(i)
        return regions

    def _find_big_segments_roi(self, roilist):
        """
        Rturns a list of 'skimage.measure._regionprops.RegionProperties' corresponding
        to the principal segments in view.

        Parameters
        ----------
        roilist : skimage.labels
            List of labels identifying all the ROIs found.

        Returns
        -------
        segments : list
            List containing the properties for the selected segments.
        """
        segments = []
        regions = measure.regionprops(roilist)
        for region in regions:
            if region.area > 10000:
                segments.append(region)
        return segments
