"""
M4 Roi utilities

This module contains utilities for the management of the ROIs in the
specific M4 case.
"""

import numpy as _np
from opticalib import typings as _ot
from skimage import measure as _meas


def single_segment_mask(image: _ot.ImageData, apply: bool = True):
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
    roilist = _meas.label(_np.invert(image.mask))
    segments = _find_big_segments_roi(roilist)
    max_area = max([seg.area for seg in segments])
    act_seg = segments[
        next(i for i, obj in enumerate(segments) if obj.area == max_area)
    ]
    maski = _np.zeros(roilist.shape, dtype=bool)
    maski[_np.where(roilist == act_seg.label)] = 1
    final_roi = _np.ma.mask_or(_np.invert(maski), image.mask)
    if apply:
        out = _np.ma.masked_array(data=image, mask=final_roi)
    else:
        out = final_roi
    return out


def adjacent_segments_mask(image: _ot.ImageData, apply: bool = True):
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
    roilist = _meas.label(_np.invert(image.mask))
    segments = _find_big_segments_roi(roilist)
    maski = mask1 = mask2 = mask3 = _np.zeros(roilist.shape, dtype=bool)
    mask1[(_np.where(roilist == segments[0].label))] = 1
    mask2[(_np.where(roilist == segments[1].label))] = 1
    mask3[(_np.where(roilist == segments[2].label))] = 1
    maski = _np.ma.mask_or(_np.ma.mask_or(mask1, mask2), mask3)
    final_roi = _np.ma.mask_or(_np.invert(maski), image.mask)
    if apply:
        out = _np.ma.masked_array(data=image, mask=final_roi)
    else:
        out = final_roi
    return out


def all_segments_mask(image: _ot.ImageData, apply: bool = True):
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
    roilist = _meas.label(_np.invert(image.mask))
    segments = _find_all_segments_roi(roilist)
    masks = [_np.zeros(roilist.shape, dtype=bool)]
    for i, seg in enumerate(segments):
        mask = _np.zeros(roilist.shape, dtype=bool)
        mask[(_np.where(roilist == seg.label))] = 1
        masks.append(mask)
    for i in range(1, len(masks)):
        masks[0] = _np.ma.mask_or(masks[0], masks[i])
    final_roi = _np.ma.mask_or(_np.invert(masks[0]), image.mask)
    if apply:
        out = _np.ma.masked_array(data=image, mask=final_roi)
    else:
        out = final_roi
    return out


def _find_all_segments_roi(roilist: _ot.ImageData):
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
    regions = _meas.regionprops(roilist)
    for i, region in enumerate(regions):
        if region.area < 10000:  # TODO - Find good criteria
            regions.pop(i)
    return regions


def _find_big_segments_roi(roilist: _ot.ImageData):
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
    regions = _meas.regionprops(roilist)
    for region in regions:
        if region.area > 10000:
            segments.append(region)
    return segments
