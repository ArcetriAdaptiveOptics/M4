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

    def automatical_roi_selection(self, image, segment_or_central_view, ref_mirror_in_or_out):
        '''
        image = np.ma.masked_array
        segment_or_central_view = 0 for segment, 1 for central
        ref_mirror_in_or_out = 0 for in, 1 for out
        '''
        if segment_or_central_view == 0:
            if ref_mirror_in_or_out == 0:
                roiList = self.roiGenerator(image)
                sx = np.zeros(np.array(roiList).shape[0])
                reference_mirror = np.zeros(np.array(roiList).shape[0])
                dx = np.zeros(np.array(roiList).shape[0])
                for i in range(np.array(roiList).shape[0]):
                    sx[i] = roiList[i][200, 100].astype(int)
                    reference_mirror[i] = roiList[i][200, 200].astype(int)
                    dx[i] = roiList[i][150, 450].astype(int)
                s = np.where(sx==0)[0][0]
                c = np.where(reference_mirror==0)[0][0]
                d = np.where(dx==0)[0][0]

                roi_sx = roiList[s]
                roi_rm = roiList[c]
                roi_dx = roiList[d]
                return roi_sx, roi_rm, roi_dx
            elif ref_mirror_in_or_out == 1:
                roiList = self.roiGenerator(image)
                sx = np.zeros(np.array(roiList).shape[0])
                central = np.zeros(np.array(roiList).shape[0])
                dx = np.zeros(np.array(roiList).shape[0])
                for i in range(np.array(roiList).shape[0]):
                    sx[i] = roiList[i][200, 100].astype(int)
                    central[i] = roiList[i][300, 300].astype(int)
                    dx[i] = roiList[i][200, 400].astype(int)
                s = np.where(sx==0)[0][0]
                c = np.where(central==0)[0][0]
                d = np.where(dx==0)[0][0]

                roi_sx = roiList[s]
                roi_central = roiList[c]
                roi_dx = roiList[d]
                return roi_sx, roi_central, roi_dx

        elif segment_or_central_view == 1:
            if ref_mirror_in_or_out == 0:
                pass
            elif ref_mirror_in_or_out == 1:
                pass

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
        