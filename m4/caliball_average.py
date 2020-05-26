'''
@author: cselmi
'''

import os
import logging
import numpy as np
from astropy.io import fits as pyfits
from m4.ground.configuration import Configuration
from m4.ground.interferometer_converter import InterferometerConverter

class Caliball():

    def __init__(self):
        """The constructor """
        self._maskthreshold = 5
        self._logger = logging.getLogger('CALIBALL:')
        self._ic = InterferometerConverter()


    @staticmethod
    def _storageFolder():
        """ Creates the path for measurement data"""
        return os.path.join(Configuration.OPD_DATA_FOLDER,
                            "Caliball")

    def createImgCubeFile(self):
        path = Caliball._storageFolder()
        total_cube = None
        for j in range(1,4):
            fold = os.path.join(path, 'test%d' %j)
            cube = self._createMeasurementCube(fold)
            if total_cube is None:
                total_cube = cube
            else:
                total_cube = np.ma.dstack((total_cube, cube))
        self._saveCube(total_cube)
        return total_cube

    def createRsImgFile(self):
        cube = self._readCube()
        mask_point = np.zeros(cube.shape[2])
        for i in range(mask_point.shape[0]):
            mask_point[i] = np.sum(cube[:,:,i].mask)
        idx = np.where(mask_point <= (min(mask_point) + self._maskthreshold))
        tot = idx[0].shape[0]
        image = np.sum(cube[:,:,idx], 3) / tot
        mask = np.prod(cube[:,:,idx].mask, 3)
        rs_img = np.ma.masked_array(image[:,:,0], mask=mask)

        fits_file_name = os.path.join(Caliball._storageFolder(), 'rs_img.fits')
        pyfits.writeto(fits_file_name, rs_img.data)
        pyfits.append(fits_file_name, rs_img.mask.astype(int))
        return rs_img

    def _createMeasurementCube(self, fold):
        cube = None
        list = os.listdir(fold)
        for i in range(len(list)-1):
            name = 'img_%04d.h5' %i
            file_name = os.path.join(fold, name)
            ima = self._ic.from4D(file_name)
            if cube is None:
                cube = ima
            else:
                cube = np.ma.dstack((cube, ima))
        return cube

    def _saveCube(self, total_cube):
        fits_file_name = os.path.join(Caliball._storageFolder(), 'Total_Cube.fits')
        pyfits.writeto(fits_file_name, total_cube.data)
        pyfits.append(fits_file_name, total_cube.mask.astype(int))

    def _readCube(self):
        file_name = os.path.join(Caliball._storageFolder(), 'Total_Cube.fits')
        hduList = pyfits.open(file_name)
        cube = np.ma.masked_array(hduList[0].data,
                                  hduList[1].data.astype(bool))
        return cube
