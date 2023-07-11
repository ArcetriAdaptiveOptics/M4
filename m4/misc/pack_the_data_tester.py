'''
Authors
  - C. Selmi: written in 2023
'''
import os
import numpy as np
from astropy.io import fits as pyfits
from m4.devices.i4d import I4D
from m4.configuration.ott_parameters import Interferometer
from m4.devices.interferometer import I4d6110

folder_test_path = '?'

def obtainDifferentTypeOfData():
    '''
    data_array = float32
    masked_ima = float64 + boolean
    masked_ima_32 = float32 + boolean
    '''
    i4d = I4D(Interferometer.i4d_IP, Interferometer.i4d_port)
    width, height, pixel_size_in_microns, data_array = i4d.takeSingleMeasurement()
    data = np.reshape(data_array, (width, height))

    interf = I4d6110()
    masked_ima = interf.acquire_phasemap()

    masked_ima_32 = masked_ima.astype(np.float32)

    return data, masked_ima, masked_ima_32


def saveData(data, masked_ima, masked_ima_32):
    fits_file_name = os.path.join(folder_test_path, 'data')
    pyfits.writeto(fits_file_name, data)

    fits_file_name = os.path.join(folder_test_path, 'masked64_maskInt')
    pyfits.writeto(fits_file_name, masked_ima.data)
    pyfits.append(fits_file_name, masked_ima.mask.astype(int))

    fits_file_name = os.path.join(folder_test_path, 'masked64_maskBool')
    pyfits.writeto(fits_file_name, masked_ima.data)
    pyfits.append(fits_file_name, masked_ima.mask)

    fits_file_name = os.path.join(folder_test_path, 'masked32_maskInt')
    pyfits.writeto(fits_file_name, masked_ima_32.data)
    pyfits.append(fits_file_name, masked_ima_32.mask.astype(int))

    fits_file_name = os.path.join(folder_test_path, 'masked32_maskBool')
    pyfits.writeto(fits_file_name, masked_ima_32.data)
    pyfits.append(fits_file_name, masked_ima_32.mask)

