"""
Authors
  - C. Selmi: written in April 2020
"""

import numpy as np


def modeRebinner(img, rebin):
    """
    Image rebinner

    Rebins a masked array image by a factor rebin.

    Parameters
    ----------
    img : masked_array
        Image to rebin.
    rebin : int
        Rebinning factor.
    
    Returns
    -------
    newImg : masked_array
        Rebinned image.
    """
    shape = img.shape
    new_shape = (shape[0]//rebin, shape[1]//rebin)
    newImg = rebin2DArray(img, new_shape)
    return newImg


def cubeRebinner(cube, rebin):
    """
    Cube rebinner

    Parameters
    ----------
    cube : ndarray
        Cube to rebin.
    rebin : int
        Rebinning factor.
    
    Returns
    -------
    newCube : ndarray
        Rebinned cube.
    """
    newCube = []
    for i in range(cube.shape[0]):
        newCube.append(modeRebinner(cube[:,:,i], rebin))
    return np.ma.dstack(newCube)


def rebin(a, *args):
    # modRB to correct for the wrong amplitude after rebin
    """rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)

    Parameters
    ----------
        a: numpy array
            vector for rebin
        *args: int
            new dimension for vector
    """
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape) / np.asarray(args)
    evList = (
        ["a.reshape("]
        + ["args[%d],factor[%d]," % (i, i) for i in range(lenShape)]
        + [")"]
        + [".sum(%d)" % (i + 1) for i in range(lenShape)]
        + ["/factor[%d]" % i for i in range(lenShape)]
    )
    print("".join(evList))
    return eval("".join(evList))


def rebin2(a, shape):
    sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
    return a.reshape(sh).mean(-1).mean(1)


# From ARTE #
def rebin2DArray(a, new_shape, sample=False):
    """
    Replacement of IDL's rebin() function for 2d arrays.
    Resizes a 2d array by averaging or repeating elements.
    New dimensions must be integral factors of original dimensions,
    otherwise a ValueError exception will be raised.
    Parameters
    ----------
    a : ndarray
        Input array.
    new_shape : 2-elements sequence
        Shape of the output array
    sample : bool
        if True, when reducing the array side elements are set
        using a nearest-neighbor algorithm instead of averaging.
        This parameter has no effect when enlarging the array.
    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array  the data are averaged,
        unless the sample parameter is set.
        If the new shape is bigger array elements are repeated.
    Raises
    ------
    ValueError
        in the following cases:
         - new_shape is not a sequence of 2 values that can be converted to int
         - new dimensions are not an integral factor of original dimensions
    NotImplementedError
         - one dimension requires an upsampling while the other requires
           a downsampling
    Examples
    --------
    >>> a = np.array([[0, 1], [2, 3]])
    >>> b = rebin(a, (4, 6)) #upsize
    >>> b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])
    >>> rebin(b, (2, 3)) #downsize
    array([[0. , 0.5, 1. ],
           [2. , 2.5, 3. ]])
    >>> rebin(b, (2, 3), sample=True) #downsize
    array([[0, 0, 1],
           [2, 2, 3]])
    """

    # unpack early to allow any 2-length type for new_shape
    m, n = map(int, new_shape)

    if a.shape == (m, n):
        return a

    M, N = a.shape

    if m <= M and n <= M:
        if (M // m != M / m) or (N // n != N / n):
            raise ValueError("Cannot downsample by non-integer factors")

    elif M <= m and M <= m:
        if (m // M != m / M) or (n // N != n / N):
            raise ValueError("Cannot upsample by non-integer factors")

    else:
        raise NotImplementedError(
            "Up- and down-sampling in different axes " "is not supported"
        )

    if sample:
        slices = [slice(0, old, float(old) / new) for old, new in zip(a.shape, (m, n))]
        idx = np.mgrid[slices].astype(int)
        return a[tuple(idx)]
    else:
        if m <= M and n <= N:
            return a.reshape((m, M // m, n, N // n)).mean(3).mean(1)
        elif m >= M and n >= M:
            return np.repeat(np.repeat(a, m / M, axis=0), n / N, axis=1)
