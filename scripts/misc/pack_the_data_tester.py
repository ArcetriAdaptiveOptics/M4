"""
Authors
  - C. Selmi: written in 2023
"""

import os
import numpy as np
from astropy.io import fits as pyfits
from m4.devices.i4d import I4D
from m4.configuration.ott_parameters import Interferometer
from m4.devices.interferometer import I4d6110
from m4.ground import read_data

folder_test_path = "/mnt/m4storage/Data/M4Data/OPTData/test/pack_the_data"


def obtainDifferentTypeOfData():
    """
    data_array = float
    masked_ima = float + boolean
    """
    i4d = I4D(Interferometer.i4d_IP, Interferometer.i4d_port)
    width, height, pixel_size_in_microns, data_array = i4d.takeSingleMeasurement()
    data = np.reshape(data_array, (width, height))

    interf = I4d6110()
    masked_ima = interf.acquire_phasemap()

    # masked_ima_32 = masked_ima.astype(float)

    return data, masked_ima


def saveData(data, masked_ima):
    fits_file_name = os.path.join(folder_test_path, "data.fits")
    pyfits.writeto(fits_file_name, data, overwrite=True)

    fits_file_name = os.path.join(folder_test_path, "masked32_maskInt.fits")
    pyfits.writeto(fits_file_name, masked_ima.data, overwrite=True)
    pyfits.append(fits_file_name, masked_ima.mask.astype(int), overwrite=True)

    fits_file_name = os.path.join(folder_test_path, "masked32_mask32.fits")
    pyfits.writeto(fits_file_name, masked_ima.data, overwrite=True)
    pyfits.append(fits_file_name, masked_ima.mask.astype(float), overwrite=True)

    fits_file_name = os.path.join(folder_test_path, "masked32_maskBool.fits")
    pyfits.writeto(fits_file_name, masked_ima.data, overwrite=True)
    pyfits.append(fits_file_name, masked_ima.mask, overwrite=True)


def openData():
    fits_file_name = os.path.join(folder_test_path, "data.fits")
    data = read_data.readFits_data(fits_file_name)

    fits_file_name = os.path.join(folder_test_path, "masked32_maskInt.fits")
    masked_ima = read_data.readFits_maskedImage(fits_file_name)

    fits_file_name = os.path.join(folder_test_path, "masked32_mask32.fits")
    masked_ima_mask32 = read_data.readFits_maskedImage(fits_file_name)

    # fits_file_name = os.path.join(folder_test_path, 'masked32_maskBool.fits')
    # masked_ima_bool = read_data.readFits_maskedImage(fits_file_name)
    return data, masked_ima, masked_ima_mask32
