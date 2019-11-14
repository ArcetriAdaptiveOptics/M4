
import numpy as np
import h5py

class InterferometerConverter(object):
    """ Class to use to convert H5 files to masked array"""

    @staticmethod
    def from4D(h5filename):
        """
        args:
            h5filename = path of h5 file to convert

        returns:
                ima = masked array image
        """
        f = h5py.File(h5filename, 'r')
        genraw = f['measurement0']['genraw']['data']
        aa = np.array(genraw)
        mask = np.zeros(aa.shape, dtype=np.bool)
        mask[np.where(aa == aa.max())] = True
        ima = np.ma.masked_array(aa * 632.8e-9, mask=mask)
        return ima
