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


def roiGenerator(ima):
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

def automatical_roi_selection(image, segment_view, ref_mirror_in):
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
    roiList = roiGenerator(image)

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

def single_segment_mask(image, apply:bool=True):
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
    segments = _find_big_segments_roi(roilist)
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

def adjacent_segments_mask(image, apply:bool=True):
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
    segments = _find_big_segments_roi(roilist)
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

def all_segments_mask(image, apply:bool=True):
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
    segments = _find_all_segments_roi(roilist)
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
    
def imgCut(img):
    """
    Cuts the image to the bounding box of the finite (non-NaN) pixels in the masked image.

    Parameters
    ----------
    image : np.ma.maskedArray
        The original masked image array.

    Returns
    -------
    cutImg = np.ma.maskedArray
        The cut image within the bounding box of finite pixels.
    """
    # Find indices of finite (non-NaN) pixels
    finite_coords = np.argwhere(np.isfinite(img))
    # If there are no finite pixels, return the original image
    if finite_coords.size == 0:
        return img
    top_left = finite_coords.min(axis=0)
    bottom_right = finite_coords.max(axis=0)
    cutImg = img[top_left[0]:bottom_right[0]+1, top_left[1]:bottom_right[1]+1] 
    return cutImg


def _find_all_segments_roi(roilist):
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

def _find_big_segments_roi(roilist):
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
