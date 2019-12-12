'''
@author: cs
'''

import logging
import numpy as np
from skimage.draw import circle
from scipy import ndimage as ndi
from m4.ground.configuration import Configuration


class ROI():
    """
    Class to be used for extracting regions of interest from the image.

    HOW TO USE IT:
    from m4.utils.roi import ROI
    roi= ROI()
    """

    def __init__(self):
        """The constructor """
        self._logger = logging.getLogger('ROI:')
        self._bigDiameter = Configuration.BIG_IMAGE_DIAMETER
        self._segmentImaDiameter = Configuration.DIAMETER_IN_PIXEL_FOR_SEGMENT_IMAGES


    def roiGenerator(self, ima):
        '''
        arg:
            ima = image (np.masked_array)

        return:
            roiList = list of the first 12 roi found in the image

        NOTA: roiList[3] = RM roi for alignement
              roiList[3] = central roi for segment

        '''
        self._logger.info('Creation of roi list')
        labels = ndi.label(np.invert(ima.mask))[0]
        #import skimage.morphology as skm
        #pro= skm.watershed(ima, markers)
        roiList = []
        for i in range(1, 13):
            maski = np.zeros(labels.shape, dtype=np.bool)
            maski[np.where(labels == i)] = 1
            final_roi = np.ma.mask_or(np.invert(maski), ima.mask)
            roiList.append(final_roi)
        return roiList

    def create_circular_mask(self, center_y, center_x, radius, imagePixels=None):
        if imagePixels is None:
            imagePixels = 512
        else:
            imagePixels = imagePixels
        mask = np.ones((imagePixels, imagePixels), dtype= bool)
        rr, cc = circle(center_y, center_x, radius)
        mask[rr,cc] = 0
        return mask

    def circularMaskForSegmentCreator(self):
        center_y = self._bigDiameter / 2
        center_x = self._bigDiameter / 2
        radius = Configuration.M4_OPTICAL_DIAMETER / 2
        big_mask = self.create_circular_mask(center_y,
                                         center_x, radius, self._bigDiameter)

        seg_center_y = np.int(self._bigDiameter/2 + Configuration.SEGMENT_DISTANCE_FROM_CENTRE)
        seg_center_x = np.int(self._bigDiameter/2)
        seg_radius = np.int(self._segmentImaDiameter / 2)
        mask = big_mask[seg_center_y - seg_radius : seg_center_y + seg_radius,
                        seg_center_x - seg_radius : seg_center_x + seg_radius]
        return mask
        