
import numpy as np
import h5py

class InterferometerConverter(object):

    @staticmethod
    def from4D(h5filename):
        f=h5py.File(h5filename, 'r')
        genraw= f['measurement0']['genraw']['data']
        aa=np.array(genraw)
        mask=np.zeros(aa.shape, dtype=np.bool)
        mask[np.where(aa==aa.max())]=True
        ima=np.ma.masked_array(aa * 632.8e-9, mask=mask)
        return ima

