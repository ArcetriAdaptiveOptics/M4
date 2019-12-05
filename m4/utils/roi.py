'''
@author: cs
'''

import logging
import numpy as np
import skimage
from scipy import ndimage as ndi


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

    def create_circular_mask(self, center_y, center_x, radius, shape):
        mask = np.ones((512, 512), dtype= bool)
        rr, cc = skimage.draw.circle(center_y, center_x, radius, shape)
        mask[rr,cc] = 0
        return mask